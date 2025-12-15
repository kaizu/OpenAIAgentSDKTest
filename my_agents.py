from datetime import datetime
from agents import Agent, function_tool


@function_tool
async def get_time() -> str:
    """現在時刻を取得する関数"""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def create_my_agent():
    return Agent(
            name="Assistant",
            model="gpt-4o-mini",
            tools=[get_time],
        )