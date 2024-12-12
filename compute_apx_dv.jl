include("src/PERC.jl")

datasets = ["00_ego-fb-combined-N.txt"   ,  "04_web_notredame.txt"     ,  "09_italian_twitter.txt"  ,"14_p2p_gnutella31.txt",
"00_hiv.txt"          ,        "05_wiki_talk.txt"   ,       "10_flickr.txt"    ,       "15_cit_hepph.txt",
"01_musae_facebook_edges.txt" , "06_web_google.txt"  ,        "11_soc_epinions.txt" ,    
"02_email_enron.txt"      ,     "07_large_twitch_edges.txt"  ,"12_soc_slashdot.txt",
"03_ca_astroph.txt"       ,     "08_web_berkstan.txt"  ,      "13_soc_pokec.txt"]

sample_size = 50000

separator = "\t"
graphs_path ="graphs/"
percolation_path = "percolation_states/"
for ds in datasets
    gf = graphs_path*ds
    perc = read_percolation(percolation_path*ds)
    n = lastindex(perc)
    ds_name = string(split(ds,".txt")[1])
    x = approx_d_max(n,perc,sample_size)
    save_d_v(ds_name,x[1],sample_size)
end