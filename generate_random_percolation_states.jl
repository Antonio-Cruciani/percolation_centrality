include("src/PERC.jl")



#datasets = ["00_hiv.txt","00_ego-fb-combined-N.txt","01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt","07_large_twitch_edges.txt","15_cit_hepph.txt","14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt","08_web_berkstan.txt"]

datasets = ["09_italian_twitter.txt","10_flickr.txt"]


directed = false

separator = "\t"

for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])
    save_percolation_array(ds_name,perc)


end