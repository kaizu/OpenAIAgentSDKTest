import logging
import streamlit as st
from dotenv import load_dotenv
from agents import Runner
from openai.types.responses import ResponseTextDeltaEvent
from my_agents import create_my_agent


# Load local environment variables for future OpenAI agent use
load_dotenv()

# Configure logging so my_agents logger outputs to console
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logging.getLogger("my_agents").setLevel(logging.INFO)

st.set_page_config(page_title="Echo Chat", page_icon="💬")
st.title("Echo Chat (準備版)")
st.caption("入力したテキストをそのまま返すシンプルなチャット。後でエージェントを組み込み予定。")

if "messages" not in st.session_state:
    st.session_state.messages = []

if "agent" not in st.session_state:
    st.session_state.agent = create_my_agent()

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
        result = Runner.run_streamed(st.session_state.agent, input=prompt, max_turns=20)
        async for event in result.stream_events():
            if event.type == "raw_response_event" and isinstance(event.data, ResponseTextDeltaEvent):
                full_response += event.data.delta
                yield event.data.delta
        st.session_state.messages.append({"role": "assistant", "content": full_response})

    with st.chat_message("assistant"):
        st.write_stream(stream_data)
