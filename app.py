import streamlit as st
from dotenv import load_dotenv
from datetime import datetime
from agents import Agent, Runner, function_tool


# Load local environment variables for future OpenAI agent use
load_dotenv()

@function_tool
async def get_time() -> str:
    """現在時刻を取得する関数"""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def decide_response(user_text: str) -> str:
    """Return a response to the user with lightweight context."""
    # Build a short history string to keep minimal context across turns.
    history = "\n".join(f"{m['role']}: {m['content']}" for m in st.session_state.messages)
    prompt = f"{history}\nuser: {user_text}" if history else user_text
    result = Runner.run_sync(st.session_state.agent, prompt)
    return result.final_output

st.set_page_config(page_title="Echo Chat", page_icon="💬")
st.title("Echo Chat (準備版)")
st.caption("入力したテキストをそのまま返すシンプルなチャット。後でエージェントを組み込み予定。")

if "messages" not in st.session_state:
    st.session_state.messages = []
    
if "agent" not in st.session_state:
    st.session_state.agent = Agent(
        name="Assistant",
        model="gpt-4o-mini",
        tools=[get_time],
    )

# Display history
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.write(message["content"])

# Chat input + echo response
if prompt := st.chat_input("メッセージを入力"):
    # User message
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.write(prompt)

    # Echo assistant message
    echo_text = decide_response(prompt)
    st.session_state.messages.append({"role": "assistant", "content": echo_text})
    with st.chat_message("assistant"):
        st.write(echo_text)
