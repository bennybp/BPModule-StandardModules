#!/usr/bin/env python3

import sys


if len(sys.argv) != 2:
   print("Usage: generate-potential-recur.py maxam")
   quit(1)



def IterateCartesian(g, am):
  if g[2] == am:
    newg = (-1, -1, -1)
  elif g[2] < (am - g[0]):
    newg = (g[0], g[1]-1, g[2]+1)
  else:
    newg = (g[0]-1, am-g[0]+1, 0)

  return newg

def NewCartesian(g, pos, n):
  l = list(g)
  l[pos] -= n
  return tuple(l)

def IsValid(g, am):
  return (g[0] >= 0 and g[0] <= am and
          g[1] >= 0 and g[1] <= am and
          g[2] >= 0 and g[2] <= am and
          g[0] + g[1] + g[2] == am)



maxl = int(sys.argv[1])

cartmap = {}

for l in range(0, maxl+1):
    g = (l, 0, 0)

    cartmap[l] = []

    while(IsValid(g, l)):
      cartmap[l].append(g)
      g = IterateCartesian(g, l)



recurmap = {}

for am,glist in cartmap.items():
    recurmap[am] = []
    for c in glist:
        newent = { "cart": c }

        #determine the direction of recurrence
        srt = sorted(c)
        for v in srt:
            if v > 0:
                break

        if v == 0:
            pos = -1
        else:
            pos = 2-list(reversed(c)).index(v)
        newent['dir'] = pos

        for d in range(0,3):
            # First and fourth terms - subtract one at position pos
            nc1 = NewCartesian(c, d, 1)
            if IsValid(nc1, am-1):
                nc1idx = cartmap[am-1].index(nc1)
            else:
                nc1idx = -1

            nc2 = NewCartesian(c, d, 2)
            if IsValid(nc2, am-2):
                nc2idx = cartmap[am-2].index(nc2)
            else:
                nc2idx = -1

            newent['idx' + str(d)] = (nc1idx, nc2idx)

        recurmap[am].append(newent)


print("RecurMap am_recur_map{")
for am,info in recurmap.items():
  print("    {{".format(am)) # opens inner vector
  for l in info:
    c = l['cart']
    print("          {{ {{ {}, {}, {} }}, {}, {{  {{ {}, {} }}, {{ {}, {} }}, {{ {}, {} }} }} }},".format(
                                                c[0], c[1], c[2],
                                                l['dir'],
                                                l['idx0'][0], l['idx0'][1],
                                                l['idx1'][0], l['idx1'][1],
                                                l['idx2'][0], l['idx2'][1]))
  print("    },") # closes inner vector pair

print("};")



 



