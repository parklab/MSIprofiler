#rus={1,2,3,4} # repeat units considered ## take this from the CONF FILE
match_score = 1 # score for a match
mismatch_score = -6 # score penalty for a mismatch
fail_score = -1 # score value to stop search at a given current_pos
min_score = 5 ##minimum score value to pass


def find_repeats_target(seq,flank_size,rus):
    cigar = seq #XX change
    bases=len(seq) # number of bases in the input sequence

    # save output as a list of lists
    out = []; exclude=set() # use sets: they are much faster with 'in'[]#np.array([],dtype='int')
    
    for ru in [rus]:
        positions_motif = range(0,ru)
        nb_positions_motif = len(positions_motif)
    # note that the flank is a range, whether python is zero-based
        not_found = True
        base = flank_size 

        while base < bases-flank_size: #and base not in exclude:
            #print base, "BASE"
            if base in exclude:
                base+=1
                continue
            elif not_found:# and base != flank_size:
               # base=test_pos+1
                test_pos=base+ru#+1flank_size
                current_pos=base
            else:
                #base=1+base+mm#test_pos#+1
                #print base
                current_pos=base
                not_found=True
                test_pos=test_pos+ru#+1#2

            pos_in_motif = 0 #update_current_pos(ru)
            score = 0; depth = 0; keep = 0
    
            max_observed_score = 0
            scores = []
            while ( (test_pos ) < (bases-flank_size) )  and  score > fail_score and test_pos not in exclude: #XX the minus one check
                #print base, "   ", score, seq[base:test_pos], depth, exclude, bases-flank_size, "RU",ru
                #print "RU current_pos  pos_in_motif  test_pos   bases"
                #print ru, current_pos, pos_in_motif, test_pos, bases #, current_pos + pos_in_motif, bases
                match = (seq[current_pos + pos_in_motif] == seq[test_pos])
        
                if match:
                    test_pos+=1
                    pos_in_motif = positions_motif[(pos_in_motif + 1) % nb_positions_motif]
                    score+=match_score
                    scores.append(score)
                    depth = 0
            
                else: # no mismatch: check for N, insertions, deletions and missense
                    score+=mismatch_score
                    scores.append(score)
                    pos_in_motif = positions_motif[(pos_in_motif + 1) % nb_positions_motif]

                    if score > fail_score and depth < 5:
                        depth+=1
                        test_pos +=1

#                        #XX if missense
#                        if cigar[test_pos] == "D":
#                            test_pos = test_pos #no avanzo
#                            #pos_in_motif =0#+=1
#                            #pos_in_motif.__next__()
#                            pos_in_motif = positions_motif[(pos_in_motif + 1) % nb_positions_motif]
#                            #run this function again ## XX see c code XXX
#                            depth+=1
#                            #partial score #XXX
#        
#                        if cigar[test_pos] == "I":
#                            test_pos +=1 # advance
#                            ## pos_in_motif DO NOT UPDATE pos_in_motif.__next__()
#                            depth+=1
    
                # keep track of the best observed score
                if score > max_observed_score:
                    max_observed_score = score
            test_pos = test_pos #- depth 
            if max_observed_score >= min_score:
                mm = scores.index(max(scores)) 
                mm = mm +ru
                out.append( [ru, base, base+mm, seq[base:base+mm+1]])#test_pos]] )
                not_found = False
                exclude.update(range(base,base+mm+1)) #np.unique(np.append(exclude,np.arange(base,test_pos)))
                test_pos = base + mm  ##X
                base = test_pos
            # increment base 
            base+=1
        else:
            base+=1
    return out
