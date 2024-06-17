include("src/PERC.jl")


datasets = [
    "00_hiv.txt",
    "00_ego-fb-combined-N.txt",
    "01_musae_facebook_edges.txt",
    "02_email_enron.txt",
    "03_ca_astroph.txt",
    "07_large_twitch_edges.txt"
    

]

directed = false

distance_metrics_stats(datasets,directed)
graph_properties(datasets,directed)
datasets = ["15_cit_hepph.txt",
    "14_p2p_gnutella31.txt",
    "11_soc_epinions.txt",
    "12_soc_slashdot.txt",
    "04_web_notredame.txt",
    "06_web_google.txt"
   ]
directed = true
distance_metrics_stats(datasets,directed)
graph_properties(datasets,directed)
