cluster.functions = makeClusterFunctionsSGE("sge-simple.tmpl") 
default.resources<-list(queue="normal.q",
                        memory = 2000,
                        pp.size = 2000)