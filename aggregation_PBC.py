import numpy as np
from collections import defaultdict

def calculate_aggregates(pdb_frame, box_length):
    def extract_coordinates(lines):
        coords = []
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    # Correct: fixed-column format
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue
        return np.array(coords)

    # Step 1: Extract each sugar’s atoms
    sugars = []
    current = []

    for line in pdb_frame:
        if line.startswith(("ATOM", "HETATM")):
            current.append(line)
        elif line.startswith(("TER", "ENDMDL")):
            coords = extract_coordinates(current)
            if len(coords) > 0:
                sugars.append(coords)
            current = []

    if current:
        coords = extract_coordinates(current)
        if len(coords) > 0:
            sugars.append(coords)

    n = len(sugars)
    contact_graph = defaultdict(set)

    # Step 2: Check pairwise distances with PBC correction
    for i in range(n):
        for j in range(i + 1, n):
            a1 = sugars[i][:, None, :]
            a2 = sugars[j][None, :, :]
            delta = a1 - a2
            delta -= box_length * np.round(delta / box_length)
            distances = np.linalg.norm(delta, axis=-1)
            if np.any(distances < 3.5):
                contact_graph[i].add(j)
                contact_graph[j].add(i)

    # Step 3: DFS to find connected components
    visited = set()
    aggregates = []

    for i in range(n):
        if i not in visited:
            stack = [i]
            group = set()
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    group.add(node)
                    stack.extend(contact_graph[node])
            aggregates.append(group)

    # Step 4: Keep aggregates with 2 or more sugars
    true_aggregates = [group for group in aggregates if len(group) >= 2]

    print("Valid aggregates (2+ sugars):", true_aggregates)

    if not true_aggregates:
        return 0, 0.0
    else:
        sizes = [len(g) for g in true_aggregates]
        return len(true_aggregates), sum(sizes) / len(sizes)
