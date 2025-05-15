from Bio import Entrez,SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class GenBankFetcher:
    def __init__(self,email,api_key):
        Entrez.email=email
        Entrez.api_key=api_key
        Entrez.tool='GenBankFetcherTool'

    def search_taxid(self,taxid):
        handle=Entrez.efetch(db="taxonomy",id=taxid,retmode="xml")
        records=Entrez.read(handle)
        org=records[0]["ScientificName"]
        print(f"Organism: {org}(TaxID: {taxid})")
        term=f"txid{taxid}[Organism]"
        handle=Entrez.esearch(db="nucleotide",term=term,usehistory="y")
        res=Entrez.read(handle)
        count=int(res["Count"])
        print(f"Found {count} records.")
        return res if count>0 else None

    def fetch_sequences(self,search_results,min_len,max_len,max_records=100):
        WebEnv=search_results["WebEnv"]
        QueryKey=search_results["QueryKey"]
        count=int(search_results["Count"])
        batch_size=500
        records=[]
        for start in range(0,min(count,max_records),batch_size):
            handle=Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=start,retmax=min(batch_size,max_records-start),webenv=WebEnv,query_key=QueryKey)
            for record in SeqIO.parse(handle,"genbank"):
                length=len(record.seq)
                if min_len<=length<=max_len:
                    records.append({
                        "accession":record.id,
                        "length":length,
                        "description":record.description
                    })
            handle.close()
        return records

def save_csv(records,filename):
    df=pd.DataFrame(records)
    df.to_csv(filename,index=False)
    print(f"CSV saved to {filename}")
    return df

def plot_lengths(df, png_file):
    df_sorted = df.sort_values("length",ascending=False)
    plt.figure(figsize=(12,6))
    plt.plot(df_sorted["accession"], df_sorted["length"],marker='o')
    plt.xticks(rotation=90,fontsize=8)
    plt.xlabel("Accession Number")
    plt.ylabel("Sequence Length")
    plt.title("GenBank Sequence Lengths")
    plt.tight_layout()
    plt.savefig(png_file)
    print(f"Plot saved to {png_file}")

def main():
    email=input("Your email: ")
    api_key=input("Your NCBI API key: ")
    taxid=input("TaxID: ")
    min_len=int(input("Minimum sequence length: "))
    max_len=int(input("Maximum sequence length: "))
    max_records=int(input("Max records to retrieve (e.g., 100): "))

    fetcher=GenBankFetcher(email, api_key)
    search=fetcher.search_taxid(taxid)
    if not search:
        print("No records found.")
        return

    records=fetcher.fetch_sequences(search, min_len, max_len, max_records)
    if not records:
        print("No records within length range.")
        return

    csv_file=f"taxid_{taxid}_report.csv"
    png_file=f"taxid_{taxid}_plot.png"
    df=save_csv(records, csv_file)
    plot_lengths(df, png_file)

if __name__=="__main__":
    main()
