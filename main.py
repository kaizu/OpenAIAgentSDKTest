from datetime import datetime
from dotenv import load_dotenv
from agents import Agent, Runner, function_tool

load_dotenv()

@function_tool
async def get_time() -> str:
    """現在時刻を取得する関数"""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def whattime():
    agent = Agent(
        name="Assistant",
        # model="gpt-4o-mini",
        model="gpt-5.2",
        tools=[get_time],
    )
    result = Runner.run_sync(agent, "いま何時？")
    print(result.final_output)

if __name__ == "__main__":
    whattime()
