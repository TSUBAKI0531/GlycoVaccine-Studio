import streamlit as st
import streamlit.components.v1 as components
from your_module import ComplexBuilder

# --- 3D 可視化コンポーネント (3Dmol.js) ---
def show_3d_viewer(pdb_text):
    """PDB データを 3D レンダリングする"""
    html_code = f"""
    <div id="container" style="height: 500px; width: 100%;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        $(function() {{
            let viewer = $3Dmol.createViewer($("#container"), {{backgroundColor: "white"}});
            viewer.addModel(`{pdb_text}`, "pdb");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.zoomTo();
            viewer.render();
        }});
    </script>
    """
    components.html(html_code, height=520)

st.set_page_config(page_title="GlycoVaccine Studio v3.0", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v3.0")

# 最初のタブ「複合体作製」のみを実装
tab1, = st.tabs(["🧬 1. 複合体作製"])

with tab1:
    st.header("Antigen-Glycan Complex Builder")
    
    # 配列および SMILES 入力
    prot_seq = st.text_area("Carrier Protein Sequence (e.g., CRM197)", height=200)
    col1, col2 = st.columns(2)
    with col1:
        linker_smi = st.text_input("Linker SMILES", "NCCCO")
    with col2:
        glycan_smi = st.text_input("Glycan SMILES", "NCCCO[C@@H]1[C@@H](NC(C)=O)[C@@H](O)O[C@H](CO)[C@H]1O")
    
    st.divider()
    
    if st.button("🛠️ 複合体を作製して表示"):
        if prot_seq:
            builder = ComplexBuilder()
            # 座標データの生成
            pdb_data = builder.build_complex_pdb(prot_seq, linker_smi, glycan_smi)
            
            st.success(f"アミノ酸数 {len(prot_seq)} のバックボーンモデルを生成しました。")
            
            # アプリ内 3D 表示
            show_3d_viewer(pdb_data)
            
            # CueMol2 用 PDB のダウンロード
            st.download_button(
                label="📥 CueMol2 用 PDB をダウンロード",
                data=pdb_data,
                file_name="antigen_complex.pdb",
                mime="chemical/x-pdb"
            )
        else:
            st.warning("キャリアタンパク質の配列を入力してください。")