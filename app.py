import streamlit as st
import pandas as pd
from your_module import ComplexBuilder, AntibodyGraftingEngine, CDRPredictor, HotSpotAnalyzer

st.set_page_config(page_title="GlycoVaccine Studio v2.0", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v2.0")

# --- Session State 管理 ---
if 'ant_seq' not in st.session_state: st.session_state.ant_seq = ""
if 'h_seq' not in st.session_state: st.session_state.h_seq = ""
if 'l_seq' not in st.session_state: st.session_state.l_seq = ""

tab1, tab2, tab3, tab4 = st.tabs(["🧬 1. 複合体作製", "🎨 2. 抗体設計", "🛡️ 3. 結合図出力", "🔥 4. Hot Spot"])

with tab1:
    st.header("🧬 Antigen-Glycan Complex Builder")
    prot = st.text_area("Carrier Protein Sequence")
    l_smi = st.text_input("Linker SMILES", "NCCCO")
    g_smi = st.text_input("Glycan SMILES", "NCCCO[C@@H]1[C@@H](NC(C)=O)[C@@H](O)O[C@H](CO)[C@H]1O")
    if st.button("Build & Export Complex"):
        builder = ComplexBuilder()
        pdb_data = builder.build_antigen_pdb(prot, l_smi, g_smi)
        st.session_state.ant_seq = prot
        st.download_button("Download PDB", pdb_data, "antigen_complex.pdb")

with tab2:
    st.header("🎨 Antibody Engineering (Trastuzumab Grafting)")
    if st.button("Graft Predicted CDRs"):
        pred = CDRPredictor()
        h_cdrs, l_cdrs = pred.predict(g_smi)
        engine = AntibodyGraftingEngine()
        st.session_state.h_seq, st.session_state.l_seq = engine.graft(h_cdrs, l_cdrs)
        st.success("CDR をトラスツズマブ FR へ移植しました。")
        st.text_area("Heavy Chain", st.session_state.h_seq)
        st.text_area("Light Chain", st.session_state.l_seq)

with tab3:
    st.header("🛡️ Output for CueMol2")
    if st.session_state.ant_seq and st.session_state.h_seq:
        builder = ComplexBuilder()
        merged = builder.merge_for_cuemol(st.session_state.ant_seq, st.session_state.h_seq, st.session_state.l_seq)
        st.download_button("Download Merged PDB for CueMol2", merged, "complex_for_cuemol.pdb")
    else:
        st.warning("Tab 1 と 2 でデータを作成してください。")

with tab4:
    st.header("🔥 Hot Spot Analysis")
    uploaded = st.file_uploader("Upload Antibody Model")
    if uploaded:
        ana = HotSpotAnalyzer()
        df = ana.analyze(uploaded.name)
        st.table(df)
        st.bar_chart(df.set_index("Residue"))