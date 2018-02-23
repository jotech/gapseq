
def seed_fix(model, fixDirection=True, rmLoops=True):
    
    mod = model.copy()
    CountDir = 0
    CountRm  = 0

    # CORRECTION IN SEED REFERENCE DATABASE
    if fixDirection:
        if "rxn05115_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn05115_c0.upper_bound=1000 # wrong direction
        if "rxn03031_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn03031_c0.lower_bound=-1000 # bidirectional http://www.rhea-db.org/reaction?id=17325
        if "rxn09272_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn09272_c0.lower_bound=0 # succinate dehydrogenase is irreversible (vmh, metacyc) 
        if "rxn10121_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn10121_c0.lower_bound=0 # Nitrate reductase Menaquinol is irreversible (vmh) 
        if "rxn00772_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn00772_c0.lower_bound=0 # ribokinase,2.7.1.15 is irreversible (vmh,metacyc)
        if "rxn01545_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn01545_c0.lower_bound=0 # Xanthosine hydrolase,3.2.2.8 is irreversible (vmh,metacyc)
        if "rxn00610_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn00610_c0.lower_bound=0 #  3.1.3.21 sn-glycerol-3-phosphate phosphohydrolase, should be irreversible (metacyc) otherwise causing loop with atpase: rxn00615 ATP glycerol 3 phosphotransferase
        if "rxn00132_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn00132_c0.lower_bound=0 # 3.1.3.5 is irreversible (metacyc,vmh) causing loop with atp rxn00134
        if "rxn00123_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn00123_c0.lower_bound=0 # 3.1.3.74 is irreversible (metacyc) causing loop with atp rxn00124
        if "rxn15029_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn15029_c0.lower_bound=0 # 3.1.3.9 is irreversible (metacyc) causing loop with atp rxn05145
        if "rxn00220_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn00220_c0.lower_bound=0 # 3.1.3.9 is irreversible (metacyc) causing loop with atp rxn05145
        if "rxn01807_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn01807_c0.lower_bound=0 # 2.7.1.35 Kinase is irreversible (metacyc,vmh)
        if "rxn01396_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn01396_c0.lower_bound=0 #  2.7.1.35 Kinase is irreversible (metacyc,vmh)
        if "rxn01671_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn01671_c0.lower_bound=0 #  irreversible (vmh)
        if "rxn00148_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn00148_c0.upper_bound=0 #  pyruvate kinase is irreversible 2.7.1.40 (metacyc,vmh)
        if "rxn00151_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn00151_c0.upper_bound=0 #  pyruvate dikinase is irreversible 2.7.9.1 (metacyc,vmh)
        if "rxn02260_0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn02260_c0.lower_bound=0 #  6.2.1.36 -- 3-hydroxypropionyl-CoA synthase is irreversible (metacyc) 
        if "rxn05054_c0" in mod.reactions:
            CountDir += 1
            mod.reactions.rxn05054_c0.lower_bound=0 #  adenosylcobinamide-phosphate synthase,6.3.1.10 irreversible (metacyc, vmh)
        #if "rxn_c0" in mod.reactions:
        #    CountDir += 1
        #    mod.reactions.rxn_c0.lower_bound=0 #  



    # CORRECTION IN SEED MODEL DATABASE
    if rmLoops:
        rmRea = set()
        if "rxn05213_c0" in mod.reactions and "rxn30467_c0" in mod.reactions: # citrate
            rmRea.add("rxn05213_c0")
        if "rxn05214_c0" in mod.reactions and "rxn30467_c0" in mod.reactions: # citrate
            rmRea.add("rxn05214_c0")
        if "rxn10816_c0" in mod.reactions and "rxn05602_c0" in mod.reactions: # lactate
            rmRea.add("rxn10816_c0")
        if "rxn10818_c0" in mod.reactions and "rxn05602_c0" in mod.reactions: # lactate
            rmRea.add("rxn10818_c0")
        if "rxn05604_c0" in mod.reactions and "rxn10818_c0" in mod.reactions: # lactate
            rmRea.add("rxn05604_c0")
        if "rxn08062_c0" in mod.reactions and "rxn05488_c0" in mod.reactions: # acetate
            rmRea.add("rxn08062_c0")
        if "rxn00145_c0" in mod.reactions and "rxn10043_c0" in mod.reactions: # cytochrome
            rmRea.add("rxn00145_c0")
        if "rxn10113_c0" in mod.reactions and "rxn10043_c0" in mod.reactions: # cytochrome
            rmRea.add("rxn10113_c0")
        if "rxn10806_c0" in mod.reactions and "rxn10043_c0" in mod.reactions: # cytochrome
            rmRea.add("rxn10806_c0")
        if "rxn05604_c0" in mod.reactions and "rxn05605_c0" in mod.reactions: # malate
            rmRea.add("rxn05604_c0")
        if "rxn10945_c0" in mod.reactions and "rxn05638_c0" in mod.reactions: # proline
            rmRea.add("rxn10945_c0")
        if "rxn05221_c0" in mod.reactions and "rxn05638_c0" in mod.reactions: # proline
            rmRea.add("rxn05221_c0")
        if "rxn05313_c0" in mod.reactions and "rxn05209_c0" in mod.reactions: # Na
            rmRea.add("rxn05313_c0")
        if "rxn05153_c0" in mod.reactions and "rxn05651_c0" in mod.reactions: # sulfate
            rmRea.add("rxn05153_c0")
        if "rxn10125_c0" in mod.reactions: # NADP transhydrogenase (direction+stoichiometry?)
            rmRea.add("rxn10125_c0")
        if "rxn05179_c0" in mod.reactions and "rxn05244_c0" in mod.reactions: # isoleucin
            rmRea.add("rxn05179_c0")
        if "rxn10138_c0" in mod.reactions and "rxn05654_c0" in mod.reactions: # succinate transport
            rmRea.add("rxn10138_c0")
        if "rxn09270_c0" in mod.reactions and "rxn10157_c0" in mod.reactions: # succinate transport, not in vmh
            rmRea.add("rxn09270_c0")
        if "rxn09270_c0" in mod.reactions and "rxn95654_c0" in mod.reactions: # succinate transport, not in vmh
            rmRea.add("rxn09270_c0")
        if "rxn00113_c0" in mod.reactions and "rxn00114_c0" in mod.reactions: # ATP carbamate phosphotransferase, not in vmh
            rmRea.add("rxn00113_c0")
        if "rxn02260_c0" in mod.reactions and "rxn02261_c0" in mod.reactions: # hydroxypropionyl coa synthase
            rmRea.add("rxn02260_c0")
        if "rxn00519_c0" in mod.reactions and "rxn00251_c0" in mod.reactions: # pep
            rmRea.add("rxn00519_c0")
        if "rxn00379_c0" in mod.reactions and "rxn09240_c0" in mod.reactions: # sulfate (H+ missing?)
            rmRea.add("rxn00379_c0")
        if "rxn00532_c0" in mod.reactions and "rxn00258_c0" in mod.reactions: # malonyl loop, not in vmh
            rmRea.add("rxn00532_c0")
        if "rxn03541_c0" in mod.reactions and "rxn05054_c0" in mod.reactions: # adenosyl cobalamide, not in vmh
            rmRea.add("rxn03541_c0")
        if "rxn16140_c0" in mod.reactions and "rxn00285_c0" in mod.reactions: # succinyl coa, not in vmh
            rmRea.add("rxn16140_c0")
        if "rxn08655_c0" in mod.reactions and "rxn00512_c0" in mod.reactions: # glycolate
            rmRea.add("rxn08655_c0")
        if "rxn08466_c0" in mod.reactions and "rxn13640_c0" in mod.reactions: # quinone+formate+H+ loop
            rmRea.add("rxn08466_c0")
        if "rxn08465_c0" in mod.reactions and "rxn13640_c0" in mod.reactions: # quinone+formate+H+ loop
            rmRea.add("rxn08465_c0")
        if "rxn10115_c0" in mod.reactions and "rxn13640_c0" in mod.reactions: # quinone+formate+H+ loop
            rmRea.add("rxn10115_c0")
        if "rxn05573_c0" in mod.reactions and "rxn05226_c0" in mod.reactions: # glc H+ transporter
            rmRea.add("rxn05773_c0")
        if "rxn05573_c0" in mod.reactions and "rxn00216_c0" in mod.reactions: # glc H+ transporter
            rmRea.add("rxn05773_c0")
        if "rxn00572_c0" in mod.reactions and "rxn10120_c0" in mod.reactions: # nitrate reductase loop, not in vmh
            rmRea.add("rxn00572_c0")
        if "rxn08907_c0" in mod.reactions and "rxn05612_c0" in mod.reactions: # melibiose transport loop, not in vmh
            rmRea.add("rxn08907_c0")
        if "rxn01134_c0" in mod.reactions: # not in vmh, causing glc-matlose-trhl loop
            rmRea.add("rxn01134_c0")
        if "rxn00134_c0" in mod.reactions and "rxn01138_c0" in mod.reactions and "rxn00131_c0" in mod.reactions : # adenine-adenosine-amp loop, not in vmh
            rmRea.add("rxn00134_c0")
        if "rxn00230_c0" in mod.reactions and "rxn00225_c0" in mod.reactions: # acetate loop, not in vmh
            rmRea.add("rxn00230_c0")
        if "rxn00606_c0" in mod.reactions and "rxn00605_c0" in mod.reactions: # trehalose loop
            rmRea.add("rxn00606_c0")
        if "rxn10150_c0" in mod.reactions and "rxn05654_c0" in mod.reactions: # succinate transporter loop
            rmRea.add("rxn10150_c0")
        if "rxn05214_c0" in mod.reactions and "rxn05213_c0" in mod.reactions: # citrate transporter loop
            rmRea.add("rxn05214_c0")
        if "rxn05740_c0" in mod.reactions and "rxn08615_c0" in mod.reactions and "rxn00695_c0" in mod.reactions : # adpglc-glycogen-glc1p loop
            rmRea.add("rxn05740_c0")
        if "rxn00339_c0" in mod.reactions and "rxn00340_c0" in mod.reactions: # 2xaspartate ammonia ligaste, not in vmh
            rmRea.add("rxn00339_c0")
        if "rxn00674_c0" in mod.reactions and "rxn00669_c0" in mod.reactions: # 2xpropionate cia ligase, not it vmh
            rmRea.add("rxn00674_c0")
        if "rxn01517_c0" in mod.reactions and "rxn01678_c0" in mod.reactions: # 2xatp phosphotransferase dUMP/dUDP
            rmRea.add("rxn01517_c0")
        if "rxn00058_c0" in mod.reactions and "rxn10043_c0" in mod.reactions: # respiratory complex 4 oxidase, reaction is not pumping protons -> loop
            rmRea.add("rxn00058_c0")
        if "rxn06106_c0" in mod.reactions and "rxn13689_c0" in mod.reactions: # respiratory complex 3, reaction is not pumping electrons -> loop
            rmRea.add("rxn06106_c0")
        if "rxn00519_c0" in mod.reactions and "rxn00515_c0" in mod.reactions: # ITP loop, not in vmh
            rmRea.add("rxn00519_c0")
        if "rxn09674_c0" in mod.reactions and "rxn05305_c0" in mod.reactions: # lysine transport loop
            rmRea.add("rxn09674_c0")
        if "rxn05573_c0" in mod.reactions and "rxn05226_c0" in mod.reactions: # glc transporter (proton symport pts loop)
            rmRea.add("rxn05573_c0")
        if "rxn05654_c0" in mod.reactions and "rxn10157_c0" in mod.reactions: # succinate transporter loop (3 and 1 proton)
            rmRea.add("rxn05654_c0")
        if "rxn00146_c0" in mod.reactions and "rxn00500_c0" in mod.reactions: # lactate cytochrome oxidureductase loop
            rmRea.add("rxn00146_c0")
        if "rxn12822_c0" in mod.reactions and "rxn05937_c0" in mod.reactions: # ferredoxin proton imbalanced, not in vmh
            rmRea.add("rxn12822_c0")
        if "rxn10860_c0" in mod.reactions and "rxn05559_c0" in mod.reactions: # formate proton symport/antiport causing loops, not in vmh
            rmRea.add("rxn10860_c0")
        if "rxn09681_c0" in mod.reactions and "rxn08535_c0" in mod.reactions: # fructose proton symport causing loops, not in vmh
            rmRea.add("rxn09681_c0")
        if "rxn00809_c0" in mod.reactions: # 1.1.1.21 this reactions is swapping nadh/nad and is proton source (+3!) could not be found in other databases, seems to be artifact causing loops
            rmRea.add("rxn00809_c0")
        if "rxn07077_c0" in mod.reactions: # wrong proton number for fad: fadh2 <-> fad + 3H ?? 
            rmRea.add("rxn07077_c0")
        if "rxn06840_c0" in mod.reactions: # 1.2.1.88  this reactions is has nad/nadh on the wrong side
            rmRea.add("rxn06840_c0")
        if "rxn_08173c0" in mod.reactions: # atpase with 4->3 proton causing imbalance
            rmRea.add("rxn08173_c0")
        if "rxn10091_c0" in mod.reactions: # cob(I)alamin adenosyltransferase: direction and compounds unclear/wrong
            rmRea.add("rxn10091_c0")
        if "rxn11268_c0" in mod.reactions: # ppi transporter from cytosol to periplasm with proton (symport) : unclear
            rmRea.add("rxn11268_c0")
        if "rxn00568_c0" in mod.reactions: # ill defined? (2) H2O[0] + (3) NAD[0] + (1) NH3[0] => (3) NADH[0] + (5) H+[0] + (1) Nitrite[0]
            rmRea.add("rxn00568_c0")
        if "rxn00569_c0" in mod.reactions: # ill defined? (2) H2O[0] + (3) NADP[0] + (1) NH      3[0] => (3) NADPH[0] + (5) H+[0] + (1) Nitrite[0]
            rmRea.add("rxn00569_c0")
        
        #if "rxn_c0" in mod.reactions and "rxn_c0" in mod.reactions: #
        #    rmRea.add("rxn_c0")

        mod.remove_reactions(rmRea) # remove unwanted reactions which could cause loops
        CountRm = len(model.reactions) - len(mod.reactions)

    print "\n\nRemoved loop causing reactions:", CountRm
    print     "Fixed direction for reactions :", CountDir, "\n"

    return(mod)
