import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.style.use('ggplot')

src="/Users/karthik/hpc_downloads/2016.09.12/outputs"

ids = [573, 1335, 562, 287, 1280, 1311]
samples = {'BSG01' : ['Klebsiella', 'Escherichia'], 'BSG03': ["Streptococcus"], 'BSG09': ["Staphylococcus"], 'BSG13': ["Streptococcus"], 'BSG16': ["Pseudomonas", "Escherichia"], 'Negctrl' :[]}

tax_df = pd.read_table("/Users/karthik/Documents/tax_dump/tax_parent.csv", sep=',', header = 0, index_col=0)
tax_df['rank'] = [i.strip() for i in tax_df['rank']]

def get_species_count(df, id, species_df, row_name, col_name):
    _t = df[df['tax_id'] == id]
    if not _t.empty:
        val = 0
        if row_name in species_df.index.values and col_name in species_df.columns.values:
            val = species_df[col_name][row_name]
        species_df = species_df.set_value(col_name, row_name, val + len(_t))

def get_counts_df(path, species_df):
    df = pd.read_table(path, sep='\t', names=["classified", "id", "tax_id", "kmers", "text"])
    name = path.split("/")
    name = name[len(name)-1]
    print(name)
    for i in ids:
        col_name = tax_df.ix[i]['name'].replace(" ","_")
        get_species_count(df, i, species_df, name, col_name)
        parent_id = tax_df.ix[i]['parent_id']
        while tax_df.ix[parent_id]['rank'] != "genus":
            parent_id = tax_df.ix[parent_id]['parent_id']
        col_name = tax_df.ix[parent_id]['name'].replace(" ","_") + "_unclassified"
        if col_name not in species_df.columns.values:
            get_species_count(df, parent_id, species_df, name, col_name)
        child = tax_df[tax_df['parent_id'] == i].index.values
        for _i in child:
            col_name = tax_df.ix[i]['name'].replace(" ","_")
            get_species_count(df, _i, species_df, name, col_name)

def format_df(df):
    for i in df.columns.values:
        if '_' not in i:
            print(i)
            
def plot_species(df):
    cols = df.columns.values
    for i in range(0, len(cols), 2):
        print(cols[i])
        print(cols[i+1])
        df[[cols[i], cols[i+1]]].plot()
        for tick in plt.gca().xaxis.get_major_ticks():
            sample_name = tick.label.get_text().split("_")[0]
            if cols[i].split("_")[0] in samples[sample_name]:                
                tick.label1.set_color("red")
                tick.label1.set_fontweight("bold")
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig("/Users/karthik/Documents/scripts/bacteremia/data/plots/species_genus/"+cols[i]+".png")
        plt.clf()
        
    
if __name__=="__main__":
    folder = "/Users/karthik/hpc_downloads/2016.09.12/outputs/"
    species_df = pd.DataFrame()
    for root, dirs, files in os.walk(folder):
        for f in files:
            get_counts_df(root+f, species_df)
    species_df = species_df.transpose()
    species_df = species_df.sort_index(axis = 1)
    species_df.index = [i.split("_")[0] for i in species_df.index.values]
    species_df.to_csv('/Users/karthik/Documents/scripts/bacteremia/data/species_plot.csv')
    species_df = species_df.drop('Undetermined_S0_L001_R1_001.kraken.full.output')
    plot_species(species_df)
