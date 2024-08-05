include("src/PERC.jl")



epsilon_list = [0.1,0.07,0.05,0.01,0.005]
sample_list = [1,2,3,4,5]


datasets = [
    "00_hiv.txt",
    "00_ego-fb-combined-N.txt",
    "01_musae_facebook_edges.txt",
    "02_email_enron.txt",
    "03_ca_astroph.txt",
    "07_large_twitch_edges.txt"
    

]
directed = false
algos = ["cmcera","era","fixed_ss"]

trials = 10

kappas = [10,25,50]

#=
j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for s in algos
        @info("ANALIZYING EPS: "*string(eps)*" ALGO "*s )
        #get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
        get_errors(ss,eps,datasets,s,trials,directed)
        get_times(ss,eps,datasets,s)
        #for kappa in kappas
        #    get_ranking_intersections(kappa,ss,eps,datasets,s,trials,directed)
        #end
    end
end

datasets = ["15_cit_hepph.txt",
    "14_p2p_gnutella31.txt",
    "11_soc_epinions.txt",
    "12_soc_slashdot.txt",
    "04_web_notredame.txt",
    "06_web_google.txt",
    ]
directed = true

j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for s in algos
        @info("ANALIZYING EPS: "*string(eps)*" ALGO "*s )
        #get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
        get_errors(ss,eps,datasets,s,trials,directed)
        get_times(ss,eps,datasets,s)
        #for kappa in kappas
        #    get_ranking_intersections(kappa,ss,eps,datasets,s,trials,directed)
        #end
    end
end


datasets = [
    "08_web_berkstan.txt"]
directed = true

j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for s in algos
        @info("ANALIZYING EPS: "*string(eps)*" ALGO "*s )
        #get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
        get_errors(ss,eps,datasets,s,trials,directed)
        get_times(ss,eps,datasets,s)
        #for kappa in kappas
        #    get_ranking_intersections(kappa,ss,eps,datasets,s,trials,directed)
        #end
    end
end
=#
# BIG
datasets = ["10_flickr.txt"]

directed = false

trials = 10

kappas = [10,25,50]


j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for s in algos
        @info("ANALIZYING EPS: "*string(eps)*" ALGO "*s )
        #get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
        #get_errors(ss,eps,datasets,s,trials,directed)
        get_times(ss,eps,datasets,s)
        #for kappa in kappas
        #    get_ranking_intersections(kappa,ss,eps,datasets,s,trials,directed)
        #end
    end
end

datasets = ["05_wiki_talk.txt","09_italian_twitter.txt","13_soc_pokec.txt"]

directed = true

j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for s in algos
        @info("ANALIZYING EPS: "*string(eps)*" ALGO "*s )
        #get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
        #get_errors(ss,eps,datasets,s,trials,directed)
        get_times(ss,eps,datasets,s)
        #for kappa in kappas
        #    get_ranking_intersections(kappa,ss,eps,datasets,s,trials,directed)
        #end
    end
end



