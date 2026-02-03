#GlycoVaccine Studio#

糖鎖抱合ワクチン設計・抗体相互作用解析の統合プラットフォーム
GlycoVaccine Studio は、腫瘍関連糖鎖（TACA）などを用いた次世代ワクチンの設計と、抗体との相互作用解析を支援するためのオールインワン・ツールです。AlphaFold 3 (AF3) による予測構造を軸に、バイオインフォマティクス的な定量的評価を提供します。

🚀 主な機能
本アプリケーションは4つの独立した解析モジュールで構成されています。

1. 🧬 抗原デザイン (Antigen Design)
キャリアタンパク質配列と糖鎖（リンカー含む）のSMILESから、AF3用の入力JSONを自動生成。

共有結合の自動特定: SMARTSパターンを用いて、タンパク質の結合部位（例：LysのNZ原子）と糖鎖の反応点を正確に紐付け。

2. 📊 モデル解析 (Model Evaluation)
AF3が出力した複数の構造モデルを一括インポート。

SASA（溶媒露出面積）解析: 糖鎖がどれだけ溶媒（抗体）に曝露しているかを算出し、ランキング。

3. 🛡️ 抗体ドッキング解析 (Paratope Analysis)
抗体-抗原複合体から結合部位（パラトープ）を自動抽出。

CDRマッピング: ANARCIを用いてIMGTナンバリングを行い、接触残基がどのCDRループ（CDR1, 2, 3）に属するかを同定。

4. 🔥 Hot Spot解析 (In silico Alanine Scanning)
結合界面の接触原子密度に基づき、結合エネルギーへの寄与が高いアミノ酸（Hot Spot）を予測。

リード抗体の最適化や変異導入試験の優先順位付けを支援。

🛠️ セットアップ
推奨環境
OS: Linux (WSL2 / Ubuntu 22.04 or later)

Python: 3.10 - 3.12

インストール手順
リポジトリのクローン

Bash
git clone https://github.com/[あなたのユーザー名]/GlycoVaccine-Studio.git
cd GlycoVaccine-Studio
仮想環境の構築とライブラリの導入

Bash
python3 -m venv venv
source venv/bin/activate
pip install streamlit pandas biopython rdkit freesasa py3Dmol
外部ツールの導入 (Optional: 抗体解析用)

Bash
sudo apt install hmmer
pip install git+https://github.com/oxpig/ANARCI.git

📖 使い方 (Quick Start)
アプリの起動

Bash
streamlit run app.py
ブラウザでアクセス: http://localhost:8501 を開きます。

ワークフローの実行:

Step 1: 「抗原デザイン」タブでJSONを作成し、AlphaFold 3 (Server等) で計算。

Step 2: 計算結果の .cif ファイルを「モデル解析」タブにドラッグ＆ドロップ。

Step 3: 最良モデルを用いて抗体とのドッキングを行い、その結果を「抗体ドッキング解析」で評価。

📂 ディレクトリ構成
Plaintext
GlycoVaccine-Studio/
├── app.py              # UI (Streamlit) フロントエンド
├── your_module.py      # 解析アルゴリズム (Backend classes)
├── README.md           # 本ファイル
└── requirements.txt    # 依存ライブラリリスト

✒️ 著者
[助田 将樹] (Ph.D. in Agriculture)

専門: 応用生命科学、免疫学、バイオインフォマティクス

背景: 抗体医薬の研究開発

📘 補足事項
本ツールは研究用プロトタイプです。AlphaFold 3の予測精度や各解析指標は、実験的な検証（表面プラズモン共鳴: SPR等）と組み合わせて解釈することを推奨します。
