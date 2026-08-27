"""Small strict Newick reader for the labelled CAPYBARA tree.

The supplied tree uses unquoted node names and numeric branch lengths. Keeping
this reader local avoids making the barcode-building pipeline depend on a large
phylogenetics package merely to traverse a named tree.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True, eq=False)
class Node:
    name: str
    length: float | None
    parent: "Node | None" = None
    children: list["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.index = 0

    def parse(self) -> Node:
        if not self.text:
            raise ValueError("Empty Newick input")
        root = self._subtree()
        self._skip_space()
        if self._peek() != ";":
            raise ValueError(f"Expected terminal ';' at character {self.index}")
        self.index += 1
        self._skip_space()
        if self.index != len(self.text):
            raise ValueError(f"Unexpected content after ';' at character {self.index}")
        return root

    def _subtree(self) -> Node:
        self._skip_space()
        if self._peek() == "(":
            self.index += 1
            children = [self._subtree()]
            while True:
                self._skip_space()
                token = self._peek()
                if token == ",":
                    self.index += 1
                    children.append(self._subtree())
                elif token == ")":
                    self.index += 1
                    break
                else:
                    raise ValueError(f"Expected ',' or ')' at character {self.index}")
            name, length = self._descriptor(required_name=False)
            node = Node(name=name, length=length, children=children)
            for child in children:
                child.parent = node
            return node
        name, length = self._descriptor(required_name=True)
        return Node(name=name, length=length)

    def _descriptor(self, required_name: bool) -> tuple[str, float | None]:
        self._skip_space()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",();":
            self.index += 1
        descriptor = self.text[start:self.index].strip()
        if not descriptor:
            if required_name:
                raise ValueError(f"Missing leaf name at character {start}")
            return "", None
        if ":" in descriptor:
            name, length_text = descriptor.rsplit(":", 1)
            name = name.strip()
            try:
                length = float(length_text)
            except ValueError as error:
                raise ValueError(f"Invalid branch length {length_text!r} at character {start}") from error
        else:
            name, length = descriptor, None
        if required_name and not name:
            raise ValueError(f"Missing leaf name at character {start}")
        return name, length

    def _skip_space(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1

    def _peek(self) -> str:
        return self.text[self.index] if self.index < len(self.text) else ""


def read_newick(path: Path) -> Node:
    return NewickParser(path.read_text(encoding="utf-8")).parse()


def preorder(root: Node):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))


def postorder(root: Node):
    stack: list[tuple[Node, bool]] = [(root, False)]
    while stack:
        node, visited = stack.pop()
        if visited:
            yield node
        else:
            stack.append((node, True))
            stack.extend((child, False) for child in reversed(node.children))

