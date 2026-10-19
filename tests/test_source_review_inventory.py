"""Attack the review inventory's claim to enumerate explicit Python callables."""

import ast
import unittest
from pathlib import Path

from scripts import assurance_audit


class SourceInventoryAttacks(unittest.TestCase):
    def test_control_flow_declarations_and_signature_lambdas_cannot_disappear(self):
        tree = ast.parse('''
if enabled:
    def hidden(default=lambda x: x):
        try:
            async def inside():
                return lambda: 1
        except Exception:
            def recovery():
                return 2
        return inside
class Container:
    if enabled:
        def method(self):
            return (lambda: (lambda: 3))
''')
        rows = list(assurance_audit.symbols(tree.body))
        actual = {(node.lineno,node.col_offset,type(node).__name__) for _,node in rows}
        expected = {(node.lineno,node.col_offset,type(node).__name__) for node in ast.walk(tree)
                    if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef,ast.Lambda))}
        self.assertEqual(actual, expected)
        self.assertEqual(len(rows), len(expected))
        names = {name for name,_ in rows}
        self.assertTrue({'hidden','hidden.inside','hidden.recovery',
                         'Container','Container.method'} <= names)
        self.assertEqual(sum('<lambda@' in name for name,_ in rows), 4)

    def test_repository_inventory_matches_independent_ast_walk_position_set(self):
        package = Path(__file__).resolve().parents[1]/'src/rocket_propulsion'
        for path in package.rglob('*.py'):
            tree = ast.parse(path.read_text(encoding='utf-8'))
            expected = {(node.lineno,node.col_offset,type(node).__name__) for node in ast.walk(tree)
                        if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef,ast.Lambda))}
            rows = list(assurance_audit.symbols(tree.body))
            actual = {(node.lineno,node.col_offset,type(node).__name__) for _,node in rows}
            with self.subTest(file=str(path)):
                self.assertEqual(actual, expected)
                self.assertEqual(len(rows), len(expected))
