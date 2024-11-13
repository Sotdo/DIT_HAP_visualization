import streamlit as st
st.set_page_config(layout="wide")

# st.logo("./images/DNA-grey-on-blue.jpg", size="large")

plot_page = st.Page("pages/plot_page.py", title="Curve plot", icon=":material/timeline:")
GOEA_page = st.Page("pages/GOEA_page.py", title="Enrichment analysis", icon=":material/search_insights:")
gene_similarity_page = st.Page("pages/gene_similarity.py", title="Gene similarity", icon=":material/difference:")

pg = st.navigation(
        {
            "Visualization": [plot_page],
            "Analysis": [GOEA_page, gene_similarity_page],
        }
    )
pg.run()