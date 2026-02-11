import streamlit as st
import streamlit.components.v1 as components
import json
from your_module import ComplexBuilder, AntibodyGraftingEngine

def show_3d_viewer(pdb_text, element_id="container"):
    """3Dmol.js ビューアー (画像ダウンロードボタン付き)"""
    pdb_json = json.dumps(pdb_text)
    html_code = f"""
    <div id="{element_id}" style="height: 500px; width: 100%; position: relative;">
        <button onclick="downloadImg('{element_id}')" style="position: absolute; top: 10px; right: 10px; z-index: 100; cursor: pointer;">🖼️ Save PNG</button>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        var viewer;
        $(function() {{
            let element = $('#{element_id}');
            viewer = $3Dmol.createViewer(element, {{backgroundColor: 'white'}});
            viewer.addModel({pdb_json}, "pdb");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.zoomTo();
            viewer.render();
        }});
        function downloadImg(id) {{
            let canvas = document.querySelector('#' + id + ' canvas');
            let link = document.createElement('a');
            link.download = 'structure.png';
            link.href = canvas.toDataURL("image/png");
            link.click();
        }}
    </script>
    """
    components.html(html_code, height=520)

st.set_page_config(page_title="GlycoVaccine Studio v4.0", layout="wide")
st.title("🧪 GlycoVaccine Studio v4.0")

# セッション管理
if 'complex_pdb' not in st.session_state: st.session_state.complex_pdb = ""
if 'target_glycan' not in st.session_state: st.session_state.target_glycan = ""

tab1, tab2 = st.tabs(["🧬 1. 複合体デザイン", "🎨 2. 抗体エンジニアリング"])

# --- Tab 1 ---
with tab1:
    st.header("Antigen-Glycan Complex Builder")
    prot = st.text_area("Carrier Protein Sequence", value="GADDVVDSSKSFVMENFSSYHGTKPGYVDSIQKGIQKPKSGTQGNYDDDWKGFYSTDNKYDAAGYSVDNENPLSGKAGGVVKVTYPGLTKVLALKVDNAETIKKELGLSLTEPLMEQVGTEEFIKRFGDGASRVVLSLPFAEGSSSVEYINNWEQAKALSVELEINFETRGKRGQDAMYEYMAQACAGNRVRR")
    col1, col2 = st.columns(2)
    with col1: l_smi = st.text_input("Linker", "NCCCO")
    with col2: g_smi = st.text_input("Glycan", "NCCCO[C@@H]1[C@@H](NC(C)=O)[C@@H](O)O[C@H](CO)[C@H]1O")

    if st.button("🛠️ 複合体を作製して表示"):
        builder = ComplexBuilder()
        st.session_state.complex_pdb = builder.build_pdb(prot, l_smi, g_smi)
        st.session_state.target_glycan = g_smi
        st.success("複合体モデルを生成しました。")

    if st.session_state.complex_pdb:
        show_3d_viewer(st.session_state.complex_pdb, "antigen_view")

# --- Tab 2 ---
with tab2:
    st.header("Antibody Grafting onto Trastuzumab")
    
    if st.session_state.target_glycan:
        if st.button("🎨 糖鎖標的に合わせて抗体を作製"):
            engine = AntibodyGraftingEngine()
            h_cdrs, l_cdrs = engine.predict_cdrs(st.session_state.target_glycan)
            h_full, l_full = engine.graft(h_cdrs, l_cdrs)
            
            st.subheader("Grafted Antibody Structure")
            # 抗体モデルの簡易 PDB 生成 (Tab 1 と同様のロジックを流用)
            builder = ComplexBuilder()
            ab_pdb = f"REMARK Grafted Antibody\n" + builder._generate_coords(h_full, offset_z=0) + "\n" + builder._generate_coords(l_full, offset_z=50)
            show_3d_viewer(ab_pdb, "antibody_view")
            
            st.download_button("📥 CueMol2 用 PDB 出力", ab_pdb, "antibody_for_cuemol.pdb")
    else:
        st.warning("先に Tab 1 で標的糖鎖の情報を入力してください。")