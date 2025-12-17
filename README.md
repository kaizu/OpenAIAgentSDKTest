## Overview
- StreamlitチャットUIでFBAエージェント・Entrezエージェントを操作するサンプルです。
- `fba_agents.py`: COBRAモデル(iJO1366)を扱うツール群（反応/代謝物/遺伝子の取得、ノックアウトなど）。
- `entrez_agent.py`: PubMed検索・要約取得ツール群（Entrez emailは環境変数から設定）。
- `app.py`: ストリーミング表示のチャットUI。

## Prerequisites
- Python 3.13+
- 必要環境変数: `OPENAI_API_KEY`（OpenAI Agents用）、`ENTREZ_EMAIL`（NCBI Entrez API用）
- 依存: `uv`ロックファイルあり。`pip install -r <(uv pip compile ...)`等、手元の手段で`pyproject.toml`/`uv.lock`に従ってインストールしてください。

## How to Run
- `.env`にキーを設定（例）
  ```
  OPENAI_API_KEY=sk-...
  ENTREZ_EMAIL=you@example.com
  ```
- インストール後、Streamlitを起動:
  ```
  streamlit run app.py
  ```
- コンソールにエージェントのツール呼び出しログが流れます。

## Notes
- Entrezツールはメールアドレス未設定時にエラーを出します。
- FBAツールは初回に`prepare_model`でモデルを読み込み直し、ノックアウト状態をリセットします。
