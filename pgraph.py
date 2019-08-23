#!/usr/bin/env python


import argparse

parser = argparse.ArgumentParser(description="Attempts to plot something sensible from a castep profile")
parser.add_argument("filename",
                    help="The profile to analyse")

parser.add_argument("-n","--network", dest="method", action="store_const",
                    const="networkx", default="plotly",
                    help="Plot a network graph rather than a Sankey diagram")

parser.add_argument("-c","--cutoff", dest="cutoff_percentage", action="store",
                    default="10",
                    help = "The percentage total time at which to disregard a subroutine")

#Parse the actual arguments
args = parser.parse_args()

#Some regexes for parsing parent and child sub names and times
import re
sub_regex   = re.compile(r"O-> (\w+)\W+\d+\W+\d+\W+(\d+\.\d+)s")
child_regex = re.compile(r"o-> (\w+)\W+\d+\W+\d+\W+(\d+\.\d+)s")

#Make a type for a profile entry
from collections import namedtuple
prof_entry = namedtuple("prof_entry",["name","time","data"])

#Parse the profile into something sensible
def read_data(filename):
    with open(filename) as fd:
        data = fd.read().split("+-----------------------------------------------------------------------------+")
    return data


translate_dict = {
    "check":"castep"
    }

if args.method == "plotly":
    #Define a tranform that we apply to all sub names
    def name_transform(name):
        n =  name.split("_")[0]
        if n in translate_dict:
            return translate_dict[n]
        else:
            return n


        
elif args.method == "networkx":
    def name_transform(name):
        return name

#Actually read the data
data = read_data(args.filename)

#Define a function to read a single entry
def parse_single_sub(data_string):
    m = sub_regex.search(data_string)
    if m:
        sub_name, sub_time = m.groups()
        #Optinally apply transformation to sub name
        sub_name = name_transform(sub_name)
        sub_time = float(sub_time)
        return prof_entry(sub_name,sub_time,data_string)
    else:
        return False
    
subs = filter(None,(parse_single_sub(d) for d in data))

#Find maximum time
total_time = max(s.time for s in subs)
min_time = total_time * float(args.cutoff_percentage) / 100.

#Filter subs based on 1% of total time
subs = filter(lambda s:s.time > min_time, subs)
sub_set = set(s.name for s in subs)

#Make a dict of node ids
node_dict = dict([s.name,i] for i,s in enumerate(subs))

label_list = [s.name for s in subs]


from collections import defaultdict
links_dict = defaultdict(float)

for s in subs:
    
    children = [prof_entry(name_transform(n),t,"") for n,t in child_regex.findall(s[2]) if name_transform(n) in sub_set]

    
    for c in children:
        if s.name != c[0]: #Skip self links
            tag = (s.name,c[0])
            links_dict[tag] += float(c[1])



source_list = []
target_list = []
value_list = []

#Finally build links lists
for s,t in links_dict:

    #Make sure link isnt too small
    if 5*links_dict[(s,t)] < min_time:
        continue
    
    if args.method == "plotly":
        source_list.append(node_dict[s])
        target_list.append(node_dict[t])
        value_list.append(links_dict[(s,t)])
    elif args.method == "networkx":
        source_list.append(s)
        target_list.append(t)
    else:
        print "Invalid method: {}".format(method)

        
if args.method == "plotly":
    import plotly.graph_objects as go
    fig = go.Figure(data=[go.Sankey(
        valueformat = ".0f",
        valuesuffix = "s",
        # Define nodes
        node = dict(
            pad = 15,
            thickness = 15,
            line = dict(color = "black", width = 0.5),
            label = label_list
        ),
        # Add links
        link = dict(
            source =  source_list,
            target =  target_list,
            value =  value_list
        ))])

    fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
    fig.show()

elif args.method == "networkx":
        
    import networkx as nx
    from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

    G = nx.DiGraph()


    for n in label_list:
        G.add_node(n)

    for s,t,v in zip(source_list,target_list,value_list):
        #G.add_edge(s,t,weight=v)
        arrow_size = max(v/1000,0.1)
        pen_width = max(v/100,0.1)
        G.add_edge(s,t,penwidth=pen_width,arrowsize=arrow_size)

    A = to_agraph(G)
    A.layout('dot')
    A.draw("{}.png".format(args.filename))
    print "Generated {}.png".format(args.filename)
