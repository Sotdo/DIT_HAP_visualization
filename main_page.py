import streamlit as st
st.set_page_config(layout="wide")
st.write("""
# DIT HAP Visualization
""")


plot_page = st.Page("pages/plot_page.py", title="plot")
test_page = st.Page("pages/data_visualization.py", title="test")

pg = st.navigation(
        {
            "Plot": [plot_page, test_page],
        }
    )
pg.run()