# -*- coding: utf-8 -*-
"""FBA Chat Streamlit app."""

import logging
import os
import streamlit as st
from dotenv import load_dotenv
from agents import Runner, Agent
from agents.tracing import set_tracing_export_api_key, set_tracing_disabled
from openai.types.responses import ResponseTextDeltaEvent
from fba_agents import create_fba_agent
from entrez_agent import create_entrez_agent


# Load local environment variables for future OpenAI agent use
load_dotenv()

# Settings in config.toml doesn't work.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logging.getLogger("fba_agents").setLevel(logging.INFO)
logging.getLogger("entrez_agents").setLevel(logging.INFO)

logger = logging.getLogger(__name__)

def create_agent():
    set_tracing_disabled(disabled=False)
    set_tracing_export_api_key(os.getenv('OPENAI_API_KEY'))
    model = os.getenv("OPENAI_MODEL", os.getenv("OPENAI_DEFAULT_MODEL"))
    fba_agent = create_fba_agent(model)
    # entrez_agent = create_entrez_agent(model)
    # orchestrator_agent = Agent(
    #     name="Orchestrator agent",
    #     instructions="""
    #     あなたは代謝工学の研究者です。文献情報を検索し、その結果をもとに反応や遺伝子をノックアウトした
    #     代謝モデルによる計算を行い、計算による予測を示します。
    #     また、モデルの情報を説明したり、計算による予測結果から仮説や次の計算条件の示唆を行います。
    #     """,
    #     tools=[
    #         fba_agent.as_tool(
    #             tool_name="simulate_with_fba",
    #             tool_description="代謝モデルの情報をユーザーにしまします。またそのモデルを使って反応や遺伝子欠損時の代謝状態の予測を行います。",
    #         ),
    #         entrez_agent.as_tool(
    #             tool_name="search_articles_with_pubmed",
    #             tool_description="PubMed検索によって文献を探索し、その要約を返します。",
    #         ),
    #     ],
    # )
    # return orchestrator_agent
    # return Agent(name="Triage agent", handoffs=[fba_agent, entrez_agent])
    return fba_agent

st.set_page_config(page_title="FBA チャットアシスタント", page_icon="💬")
st.title("FBA チャットアシスタント")
st.caption("代謝モデルに質問し、ツール実行をストリーミングで確認できます。")

if "messages" not in st.session_state:
    st.session_state.messages = []

if "agent" not in st.session_state:
    # st.session_state.agent = create_fba_agent()
    # st.session_state.agent = create_entrez_agent()
    st.session_state.agent = create_agent()

# Display history
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.write(message["content"])

# Chat input + echo response
if user_input := st.chat_input("メッセージを入力"):
    # User message
    st.session_state.messages.append({"role": "user", "content": user_input})
    with st.chat_message("user"):
        st.write(user_input)

    async def stream_data():
        full_response = ""
        history = "\n".join(f"{m['role']}: {m['content']}" for m in st.session_state.messages)
        prompt = f"{history}\nuser: {user_input}" if history else user_input
        result = Runner.run_streamed(st.session_state.agent, input=prompt, max_turns=30)
        async for event in result.stream_events():
            if event.type == "raw_response_event":
                if isinstance(event.data, ResponseTextDeltaEvent):
                    full_response += event.data.delta
                    yield event.data.delta

        st.session_state.messages.append({"role": "assistant", "content": full_response})

    with st.chat_message("assistant"):
        st.write_stream(stream_data)
