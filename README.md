## Overview

- StreamlitチャットUIでFBAエージェント・Entrezエージェントを操作するサンプルです。

## Prerequisites

- Python 3.13+

- 必要環境変数: `OPENAI_API_KEY`（OpenAI Agents用）、`ENTREZ_EMAIL`（NCBI Entrez API用）

## How to Run

- `.env`にキーを設定（例）

  ```
  OPENAI_API_KEY=sk-...
  ENTREZ_EMAIL=you@example.com
  OPENAI_MODEL=gpt-5.2
  ```

- インストール後、Streamlitを起動:

  ```
  uv run streamlit run app.py
  ```

- `http://localhost:8501/`にアクセス