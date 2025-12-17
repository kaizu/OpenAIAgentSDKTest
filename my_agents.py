import logging
import re
from typing import TypedDict, Any
from agents import Agent, function_tool, RunContextWrapper
import cobra

logger = logging.getLogger(__name__)


def my_custom_error_function(context: RunContextWrapper[Any], error: Exception) -> str:
    """A custom function to provide a user-friendly error message."""
    logger.info(f"my_custom_error_function called with error: {error}")
    print(f"A tool call failed with the following error: {error}")
    return f"An error occurred while running the tool. Please try again. Error: {str(error)}"

model: cobra.core.model.Model | None = None
solution: cobra.core.solution.Solution | None = None

@function_tool
async def prepare_model() -> bool:
    """代謝モデルを準備し、計算結果を取得できる状態にします。準備に成功したかどうかを返します。"""
    logger.info("prepare_model called")
    global model
    model = cobra.io.load_model('iJO1366')
    global solution
    solution = model.optimize()
    return True

@function_tool
async def get_objective_value() -> float:
    """代謝モデルの計算結果から目的関数の値を返します。目的関数の値は代謝モデルにおける増殖速度を意味します。モデルが準備できていない場合はエラーとなります。"""
    logger.info("get_objective_value called")
    if solution is None:
        raise RuntimeError("モデルの計算結果が準備できていません。")
    return solution.objective_value

@function_tool
async def get_flux(reaction_id: str) -> float:
    """代謝モデルの計算結果から反応の流束、すなわち速度を返します。
    モデルが準備できていない場合はエラーとなります。
    
    Args:
        reaction_id: 流束を取得したい反応の名前。
    """
    logger.info(f"get_flux called for reaction_id={reaction_id}")
    if solution is None:
        raise RuntimeError("モデルの計算結果が準備できていません。")
    elif reaction_id not in solution.fluxes:
        raise ValueError(f"{reaction_id}という反応は存在しません。")
    return solution.fluxes[reaction_id]

@function_tool
async def get_reaction_ids(pattern: str = ".*") -> list[str]:
    """代謝モデルに含まれる反応IDのリストを返します。モデルが準備できていない場合はエラーとなります。
    
    Args:
        pattern: 正規表現。これにマッチするIDのみ返します（デフォルトは全件）。
    """
    logger.info(f"get_reaction_ids called with pattern={pattern}")
    if model is None:
        raise RuntimeError("モデルが準備できていません。")
    try:
        regex = re.compile(pattern)
    except re.error as e:
        raise ValueError(f"正規表現が不正です: {e}") from e
    return sorted([r.id for r in model.reactions if regex.search(r.id)])

@function_tool
async def get_metabolite_ids(pattern: str = ".*") -> list[str]:
    """代謝モデルに含まれる代謝物IDのリストを返します。モデルが準備できていない場合はエラーとなります。
    
    Args:
        pattern: 正規表現。これにマッチするIDのみ返します（デフォルトは全件）。
    """
    logger.info(f"get_metabolite_ids called with pattern={pattern}")
    if model is None:
        raise RuntimeError("モデルが準備できていません。")
    try:
        regex = re.compile(pattern)
    except re.error as e:
        raise ValueError(f"正規表現が不正です: {e}") from e
    return sorted([m.id for m in model.metabolites if regex.search(m.id)])

@function_tool
async def get_gene_ids(pattern: str = ".*") -> list[str]:
    """代謝モデルに含まれる遺伝子IDのリストを返します。モデルが準備できていない場合はエラーとなります。
    
    Args:
        pattern: 正規表現。これにマッチするIDのみ返します（デフォルトは全件）。
    """
    logger.info(f"get_gene_ids called with pattern={pattern}")
    if model is None:
        raise RuntimeError("モデルが準備できていません。")
    try:
        regex = re.compile(pattern)
    except re.error as e:
        raise ValueError(f"正規表現が不正です: {e}") from e
    return sorted([g.id for g in model.genes if regex.search(g.id)])

class ReactionInfo(TypedDict):
    id: str
    name: str
    formula: str
    lower_bound: float
    upper_bound: float

class MetaboliteInfo(TypedDict):
    id: str
    name: str
    formula: str

# class GeneInfo(TypedDict):
#     id: str
#     name: str

@function_tool(failure_error_function=my_custom_error_function)
async def get_reaction_info(reaction_id: str) -> ReactionInfo:
    """代謝モデルに含まれる反応のIDから反応の情報を返します。情報はid, name, formula, lower_bound, upper_boundからなり、それぞれID、より詳細な名前、反応式、反応が取り得る流束の最小値と最大値を表します。
    モデルが準備できていない場合はエラーとなります。
    
    Args:
        reaction_id: 情報を取得したい反応のID。    
    """
    logger.info(f"get_reaction_info called for reaction_id={reaction_id}")
    if model is None:
        raise RuntimeError("モデルの準備できていません。")
    reaction = model.reactions.get_by_id(reaction_id)
    return ReactionInfo(id=reaction.id, name=reaction.name, formula=reaction.reaction, lower_bound=reaction.lower_bound, upper_bound=reaction.upper_bound)

@function_tool(failure_error_function=my_custom_error_function)
async def get_metabolite_info(metabolite_id: str) -> MetaboliteInfo:
    """代謝モデルに含まれる代謝物のIDから情報を返します。情報はid, name, formulaからなり、それぞれID、より詳細な名前、化学式を表します。
    モデルが準備できていない場合はエラーとなります。
    
    Args:
        metabolite_id: 情報を取得したい代謝物のID。
    """
    logger.info(f"get_metabolite_info called for metabolite_id={metabolite_id}")
    if model is None:
        raise RuntimeError("モデルの準備できていません。")
    metabolite = model.metabolites.get_by_id(metabolite_id)
    return MetaboliteInfo(id=metabolite.id, name=metabolite.name, formula=metabolite.formula)

# @function_tool(failure_error_function=my_custom_error_function)
# async def get_gene_info(gene_id: str) -> GeneInfo:
#     """代謝モデルに含まれる遺伝子IDから情報を返します。情報はid, nameからなり、それぞれIDとより詳細な名前を表します。
#     モデルが準備できていない場合はエラーとなります。
    
#     Args:
#         gene_id: 情報を取得したい遺伝子のID。
#     """
#     logger.info(f"get_gene_info called for gene_id={gene_id}")
#     if model is None:
#         raise RuntimeError("モデルの準備できていません。")
#     gene = model.genes.get_by_id(gene_id)
#     return GeneInfo(id=gene.id, name=gene.name)

@function_tool(failure_error_function=my_custom_error_function)
async def get_metabolite_associated_reaction_ids(metabolite_id: str) -> list[str]:
    """代謝モデルに含まれる代謝物のIDからそれを基質、もしくは生成物とする反応のIDのリストを返します。
    モデルが準備できていない場合はエラーとなります。
    
    Args:
        metabolite_id: 関連する反応を取得したい代謝物のID。
    
    Returns:
        reaction_ids: 与えられた代謝物を基質か生成物とする反応のIDのリスト。
    """
    logger.info(f"get_metabolite_associated_reactions called for reaction_id={metabolite_id}")
    if model is None:
        raise RuntimeError("モデルの準備できていません。")
    metabolite = model.metabolites.get_by_id(metabolite_id)
    return sorted([r.id for r in metabolite.reactions])

@function_tool(failure_error_function=my_custom_error_function)
async def get_reaction_associated_gene_ids(reaction_id: str) -> list[str]:
    """与えられたIDの反応を制御する酵素の遺伝子IDのリストを返します。モデルが準備できていない場合はエラーとなります。
    
    Args:
        reaction_id: 関連遺伝子を取得したい反応のID。
    """
    logger.info(f"get_reaction_associated_genes called for reaction_id={reaction_id}")
    if model is None:
        raise RuntimeError("モデルの準備できていません。")
    reaction = model.reactions.get_by_id(reaction_id)
    return sorted([gene_id for gene_id in reaction.gpr._genes])

@function_tool(failure_error_function=my_custom_error_function)
async def get_gene_associated_reaction_ids(gene_id: str) -> list[str]:
    """与えられた酵素遺伝子IDに関連する反応のIDリストを返します。モデルが準備できていない場合はエラーとなります。
    
    Args:
        gene_id: 関連する反応を取得したい遺伝子のID。
    """
    logger.info(f"get_gene_associated_reaction_ids called for gene_id={gene_id}")
    if model is None:
        raise RuntimeError("モデルの準備できていません。")
    gene = model.genes.get_by_id(gene_id)
    return sorted([r.id for r in gene.reactions])

def create_my_agent():
    logger.info("create_my_agent called")
    return Agent(
            name="Assistant",
            # model="gpt-4o-mini",
            model="gpt-5.2",
            tools=[
                get_objective_value,
                prepare_model,
                get_flux,
                get_reaction_ids,
                get_metabolite_ids,
                get_gene_ids,
                get_reaction_info,
                get_metabolite_info,
                get_gene_info,
                get_metabolite_associated_reaction_ids,
                get_reaction_associated_gene_ids,
                get_gene_associated_reaction_ids,
                ],
            instructions="""
            あなたは代謝シミュレーションの専門家です。ステップ毎に行ったことと結果を簡潔に説明して下さい。
            モデルとその結果に関する質問に対して、ツールを使って取得した情報に忠実に答えること。
            モデルから得た文字列情報については'\text'などの変更を加えずに可能な限りそのまま返すこと。
            """,
        )
