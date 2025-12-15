from datetime import datetime

import streamlit as st
from dotenv import load_dotenv


# Load local environment variables for future OpenAI agent use
load_dotenv()

st.set_page_config(page_title="Echo Chat", page_icon="💬")
st.title("Echo Chat (準備版)")
st.caption("入力したテキストをそのまま返すシンプルなチャット。後でエージェントを組み込み予定。")

if "messages" not in st.session_state:
    st.session_state.messages = []

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
    echo_text = prompt
    st.session_state.messages.append({"role": "assistant", "content": echo_text})
    with st.chat_message("assistant"):
        st.write(echo_text)
