import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os
from custom_upgma import parse_fasta, calculate_distance_matrix, run_upgma, draw_dendrogram

st.set_page_config(
    page_title="Phylogenetic Tree Builder",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .reportview-container {
        background: #0e1117;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: transparent;
        border-radius: 4px 4px 0px 0px;
        gap: 1px;
        padding-top: 10px;
        padding-bottom: 10px;
        font-weight: 600;
        font-size: 16px;
    }
    .stTabs [aria-selected="true"] {
        border-bottom: 3px solid #ff4b4b;
        color: #ff4b4b;
    }
</style>
""", unsafe_allow_html=True)

PRESET_SEQUENCES = {
}

if "custom_sequences" not in st.session_state:
    st.session_state.custom_sequences = {}

st.title("Interactive UPGMA Phylogenetic Tree Builder")
st.markdown("Compare species DNA sequences, adjust penalties, and visualize evolutionary trees instantly.")

st.sidebar.header("Comparison Bag & Controls")

uploaded_file = st.sidebar.file_uploader(
    "Upload custom FASTA (.fasta) file:",
    type=["fasta", "fa"]
)

if uploaded_file is not None:
    temp_path = "temp_upload.fasta"
    with open(temp_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    try:
        parsed_seqs = parse_fasta(temp_path)
        if parsed_seqs:
            st.session_state.custom_sequences.update(parsed_seqs)
            st.sidebar.success(f"Successfully loaded {len(parsed_seqs)} custom sequences!")
        else:
            st.sidebar.error("No valid sequences found in the file.")
    except Exception as e:
        st.sidebar.error(f"Error parsing file: {e}")
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)

if st.session_state.custom_sequences:
    if st.sidebar.button("Clear Custom Sequences"):
        st.session_state.custom_sequences = {}
        st.rerun()

full_bag = {**PRESET_SEQUENCES, **st.session_state.custom_sequences}

selected_species = st.sidebar.multiselect(
    "Pick species to compare:",
    options=list(full_bag.keys()),
    default=list(PRESET_SEQUENCES.keys())[:4]
)

st.sidebar.markdown("---")
st.sidebar.subheader("Algorithm Settings")

gap_penalty = st.sidebar.slider(
    "Gap (Insertion/Deletion) Penalty",
    min_value=0.1,
    max_value=3.0,
    value=1.0,
    step=0.1,
    help="Cost applied for inserting or deleting a base."
)

mismatch_penalty = st.sidebar.slider(
    "Mismatch Penalty",
    min_value=0.1,
    max_value=3.0,
    value=1.0,
    step=0.1,
    help="Cost applied for mismatching bases."
)

if len(selected_species) < 2:
    st.info("Please select at least 2 species from the Comparison Bag in the sidebar to build a tree!")
else:
    active_sequences = {sp: full_bag[sp] for sp in selected_species}
    matrix = calculate_distance_matrix(active_sequences, gap_penalty=gap_penalty, mismatch_penalty=mismatch_penalty)
    root_node = run_upgma(matrix)
    tab1, tab2, tab3 = st.tabs([
        "Selected Sequences", 
        "Distance Matrix ($D_0$)", 
        "Evolutionary Tree (Dendrogram)"
    ])
    
    with tab1:
        st.subheader("Selected Sequences in comparison bag")
        seq_data = [{"Species": sp, "Length (bp)": len(seq), "Sequence": seq} for sp, seq in active_sequences.items()]
        df_seqs = pd.DataFrame(seq_data)
        st.dataframe(df_seqs, use_container_width=True, hide_index=True)
        
        st.subheader("Raw FastA Visualizations")
        for sp, seq in active_sequences.items():
            st.markdown(f"**`>{sp}`**")
            st.code(seq, language="dna")
            
    with tab2:
        st.subheader("Calculated Distance Matrix ($D_0$)")
        st.markdown("Shows the normalized Levenshtein distance between each pair. `0.0` represents a perfect match.")
        
        df_matrix = pd.DataFrame(matrix)
        
        styled_df = df_matrix.style.background_gradient(cmap="viridis_r", axis=None).format("{:.3f}")
        st.dataframe(styled_df, use_container_width=True)
        
    with tab3:
        st.subheader("UPGMA Phylogenetic Tree Dendrogram")
        st.markdown("The tree is constructed using Unweighted Pair Group Method with Arithmetic Mean (UPGMA).")
        
        fig = draw_dendrogram(root_node, show=False)
        
        st.pyplot(fig)
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Total Tree Height (Max Distance)", f"{root_node.height:.3f}")
        with col2:
            st.metric("Root Node Cluster", root_node.name)
