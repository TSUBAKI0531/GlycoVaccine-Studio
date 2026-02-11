import streamlit as st
import streamlit.components.v1 as components
from your_module import ComplexBuilder

def show_3d_viewer(pdb_text):
    """3Dmol.js で Cartoon (Ribbon) 表示を行う"""
    # JS 内で PDB データを安全に扱うために改行を処理
    pdb_escaped = pdb_text.replace("\n", "\\n")
    html_code = f"""
    <div id="container" style="height: 500px; width: 100%;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        $(function() {{
            let viewer = $3Dmol.createViewer($("#container"), {{backgroundColor: "white"}});
            let data = `{pdb_escaped}`;
            viewer.addModel(data, "pdb");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.zoomTo();
            viewer.render();
        }});
    </script>
    """
    components.html(html_code, height=520)

st.set_page_config(page_title="GlycoVaccine Studio", layout="wide")
st.title("🧪 GlycoVaccine Studio v3.1")

tab1, = st.tabs(["🧬 1. 複合体作製"])

with tab1:
    st.header("Antigen-Glycan Complex Builder")
    prot_seq = st.text_area("Carrier Protein Sequence", height=150)
    col1, col2 = st.columns(2)
    with col1:
        l_smi = st.text_input("Linker SMILES", "NCCCO")
    with col2:
        g_smi = st.text_input("Glycan SMILES", "NCCCO[C@@H]1[C@@H](NC(C)=O)[C@@H](O)O[C@H](CO)[C@H]1O")
    
    if st.button("🛠️ 複合体を作製して表示"):
        if prot_seq:
            builder = ComplexBuilder()
            pdb_data = builder.build_complex_pdb(prot_seq, l_smi, g_smi)
            st.success("螺旋構造モデルを生成しました。")
            
            # アプリ内での表示
            show_3d_viewer(pdb_data)
            
            # ダウンロード
            st.download_button("📥 CueMol2 用 PDB をダウンロード", pdb_data, "antigen_complex.pdb")