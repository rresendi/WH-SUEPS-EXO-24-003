# making sure refhlt is passed in both mc and data
if passDen:
    if passDen and fillvar is not None:
        histos[var + "_den"].Fill(fillvar)
        if passHLT and lepton_matched:
            if getattr(ev, refhlt, False):
                histos[var + "_num"].Fill(fillvar)

# only applying refhlt to data          
if passDen:                                                                                                                                                               
    # Include refHLT in every denominator when data == "data"                                                                                                             
    if passDen and fillvar is not None:                                                                                                                                   
        histos[var + "_den"].Fill(fillvar)                                                                                                                                
        if passHLT and lepton_matched:                                                                                                                                     
            if data == "data":                                                                                                                                             
                if getattr(ev, refhlt, False):                                                                                                                             
                    histos[var + "_num"].Fill(fillvar)                                                                                                                     
                else:                                                                                                                                                          
                    histos[var + "_num"].Fill(fillvar)          
