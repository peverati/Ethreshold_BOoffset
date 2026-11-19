"""Key-by-key comparison of two L = 7 term-list pickles (written 2026-11-16).  Usage: python check_L7_terms.py OLD NEW"""
import pickle, sys
old = pickle.load(open(sys.argv[1], "rb"))
new = pickle.load(open(sys.argv[2], "rb"))
print(f"old {len(old)} keys, new {len(new)} keys; only old {len(set(old)-set(new))}, only new {len(set(new)-set(old))}")
diff = [k for k in set(old) & set(new) if repr(old[k]) != repr(new[k])]
print(f"{len(diff)} common keys with different term lists")
for k in diff[:10]: print("  ", k[0], "".join(k[1]), "".join(k[2]), "len old/new", len(old[k]) if hasattr(old[k], '__len__') else '-', len(new[k]) if hasattr(new[k], '__len__') else '-')
