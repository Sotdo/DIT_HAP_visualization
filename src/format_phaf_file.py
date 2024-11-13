# %%
import pandas as pd
from goatools.obo_parser import GODag
# %%
phaf_dag = GODag("/data/c/yangyusheng_optimized/DIT_HAP_visualization/Streamlit_DEseq2/references/pombase_annotation/20241001/fypo-simple.obo", optional_attrs=["defn", "relationship"])
phaf = pd.read_csv("/data/c/yangyusheng_optimized/DIT_HAP_visualization/Streamlit_DEseq2/references/pombase_annotation/20241001/pombase-2024-10-01.phaf", sep="\t").query("(`Allele type` == 'deletion' or `Allele type` == 'disruption') and Condition.str.contains('FYECO:0000005')")

# %%
def assign_FYPO_name(FYPO_ID, FYPO_dag):
    if FYPO_ID in FYPO_dag:
        return phaf_dag[FYPO_ID].name
    else:
        return "No record in simple FYPO"

phaf["DB"] = "PomBase"
phaf["DB_Object_ID"] = phaf["Gene systematic ID"]
phaf["DB_Object_Symbol"] = phaf["Gene symbol"]
phaf["Qualifier"] = ""
phaf["GO_ID"] = phaf["FYPO ID"]
phaf["DB:Reference"] = phaf["Reference"]
phaf["Evidence"] = phaf["Evidence"]
phaf["With"] = ""
phaf["Aspect"] = "FYPO"
phaf["DB_Object_Name"] = phaf["FYPO ID"].apply(assign_FYPO_name, FYPO_dag=phaf_dag)
phaf["Synonym"] = ""
phaf["DB_Object_Type"] = "protein"
phaf["Taxon"] = "taxon:4896"
phaf["Date"] = phaf["Date"].str.replace("-", "")
phaf["Assigned_By"] = phaf["#Database name"]
phaf["Annotation_Extension"] = phaf["Extension"]
phaf["Gene_Product_Form_ID"] = ""

reformat_phaf = phaf[["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence", "With", "Aspect",
             "DB_Object_Name", "Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"]].copy()

# %%
with open("../references/pombase_annotation//phaf_go_style_gaf.tsv", "w") as f:
    f.write("!gaf-version: 2.2\n!generated-by: PomBase\n!date-generated: 2024-09-30T21:55\n!URL: https://www.pombase.org\n!contact: helpdesk@pombase.org\n")
    
reformat_phaf.to_csv("../references/pombase_annotation/phaf_go_style_gaf.tsv", sep="\t", index=False, header=False, mode="a")
# %%
