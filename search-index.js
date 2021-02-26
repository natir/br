var searchIndex = JSON.parse('{\
"br":{"doc":"","i":[[0,"cli","br","",null,null],[5,"i82level","br::cli","",null,[[],[["option",4],["level",4]]]],[0,"correct","br","",null,null],[0,"exist","br::correct","",null,null],[0,"one","br::correct::exist","",null,null],[4,"ScenarioOne","br::correct::exist::one","",null,null],[13,"I","","",0,null],[13,"S","","",0,null],[13,"D","","",0,null],[3,"ScenarioOneIter","","",null,null],[6,"One","","",null,null],[0,"two","br::correct::exist","",null,null],[4,"ScenarioTwo","br::correct::exist::two","",null,null],[13,"II","","",1,null],[13,"IS","","",1,null],[13,"SS","","",1,null],[13,"SD","","",1,null],[13,"DD","","",1,null],[13,"ICI","","",1,null],[13,"ICS","","",1,null],[13,"ICD","","",1,null],[13,"SCI","","",1,null],[13,"SCS","","",1,null],[13,"SCD","","",1,null],[13,"DCI","","",1,null],[13,"DCD","","",1,null],[3,"ScenarioTwoIter","","",null,null],[6,"Two","","",null,null],[8,"Scenario","br::correct::exist","",null,null],[10,"init","","",2,[[]]],[10,"c","","",2,[[]]],[10,"apply","","",2,[[["boxkmerset",6]],["option",4]]],[10,"correct","","",2,[[["boxkmerset",6]]]],[11,"get_score","","",2,[[["boxkmerset",6]]]],[11,"one_more","","",2,[[["boxkmerset",6]]]],[3,"Exist","","",null,null],[11,"new","","",3,[[["boxkmerset",6]]]],[0,"gap_size","br::correct","",null,null],[3,"GapSize","br::correct::gap_size","",null,null],[11,"new","","",4,[[["boxkmerset",6]]]],[11,"ins_sub_correction","","",4,[[],["option",4]]],[0,"graph","br::correct","",null,null],[3,"Graph","br::correct::graph","",null,null],[11,"new","","",5,[[["boxkmerset",6]]]],[0,"greedy","br::correct","",null,null],[3,"Greedy","br::correct::greedy","",null,null],[11,"new","","",6,[[["boxkmerset",6]]]],[8,"Corrector","br::correct","",null,null],[10,"valid_kmer","","",7,[[],["boxkmerset",6]]],[10,"correct_error","","",7,[[],["option",4]]],[11,"k","","",7,[[]]],[11,"correct","","",7,[[],["vec",3]]],[0,"error","br","",null,null],[4,"Error","br::error","All error produce by Pcon",null,null],[13,"Cli","","See enum [Cli]",8,null],[13,"IO","","See enum [IO]",8,null],[13,"CantComputeAbundance","","",8,null],[4,"Cli","","Error emmit durring Cli parsing",null,null],[13,"NotSameNumberOfInAndOut","","Number of inputs and outputs must be the same",9,null],[13,"NoSolidityNoKmer","","",9,null],[13,"KmerSolidNeedK","","",9,null],[4,"IO","","Error emmit when pcon try to work with file",null,null],[13,"CantCreateFile","","We can\'t create file. In C binding it\'s equal to 0",10,null],[13,"CantOpenFile","","We can\'t open file. In C binding it\'s equal to 1",10,null],[13,"ErrorDurringWrite","","Error durring write in file. In C binding it\'s equal to 2",10,null],[13,"ErrorDurringRead","","Error durring read file. In C binding it\'s equal to 3",10,null],[13,"NoError","","No error, this exist only for C binding it\'s the value of …",10,null],[0,"set","br","",null,null],[0,"hash","br::set","",null,null],[3,"Hash","br::set::hash","",null,null],[11,"new","","",11,[[["vec",3]]]],[0,"pcon","br::set","",null,null],[3,"Pcon","br::set::pcon","",null,null],[11,"new","","",12,[[["solid",3]]]],[8,"KmerSet","br::set","",null,null],[10,"get","","",13,[[]]],[10,"k","","",13,[[]]],[6,"BoxKmerSet","","",null,null],[5,"run_correction","br","",null,[[["box",3],["vec",3]],["result",6]]],[5,"build_methods","","",null,[[["boxkmerset",6],["option",4],["vec",3]],[["box",3],["vec",3]]]],[5,"set_nb_threads","","Set the number of threads use by count step",null,[[]]],[11,"from","br::correct::exist::one","",0,[[]]],[11,"into","","",0,[[]]],[11,"to_owned","","",0,[[]]],[11,"clone_into","","",0,[[]]],[11,"borrow","","",0,[[]]],[11,"borrow_mut","","",0,[[]]],[11,"try_from","","",0,[[],["result",4]]],[11,"try_into","","",0,[[],["result",4]]],[11,"type_id","","",0,[[],["typeid",3]]],[11,"init","","",0,[[]]],[11,"deref","","",0,[[]]],[11,"deref_mut","","",0,[[]]],[11,"drop","","",0,[[]]],[11,"vzip","","",0,[[]]],[11,"from","","",14,[[]]],[11,"into","","",14,[[]]],[11,"into_iter","","",14,[[]]],[11,"to_owned","","",14,[[]]],[11,"clone_into","","",14,[[]]],[11,"borrow","","",14,[[]]],[11,"borrow_mut","","",14,[[]]],[11,"try_from","","",14,[[],["result",4]]],[11,"try_into","","",14,[[],["result",4]]],[11,"type_id","","",14,[[],["typeid",3]]],[11,"init","","",14,[[]]],[11,"deref","","",14,[[]]],[11,"deref_mut","","",14,[[]]],[11,"drop","","",14,[[]]],[11,"vzip","","",14,[[]]],[11,"from","br::correct::exist::two","",1,[[]]],[11,"into","","",1,[[]]],[11,"to_owned","","",1,[[]]],[11,"clone_into","","",1,[[]]],[11,"borrow","","",1,[[]]],[11,"borrow_mut","","",1,[[]]],[11,"try_from","","",1,[[],["result",4]]],[11,"try_into","","",1,[[],["result",4]]],[11,"type_id","","",1,[[],["typeid",3]]],[11,"init","","",1,[[]]],[11,"deref","","",1,[[]]],[11,"deref_mut","","",1,[[]]],[11,"drop","","",1,[[]]],[11,"vzip","","",1,[[]]],[11,"from","","",15,[[]]],[11,"into","","",15,[[]]],[11,"into_iter","","",15,[[]]],[11,"to_owned","","",15,[[]]],[11,"clone_into","","",15,[[]]],[11,"borrow","","",15,[[]]],[11,"borrow_mut","","",15,[[]]],[11,"try_from","","",15,[[],["result",4]]],[11,"try_into","","",15,[[],["result",4]]],[11,"type_id","","",15,[[],["typeid",3]]],[11,"init","","",15,[[]]],[11,"deref","","",15,[[]]],[11,"deref_mut","","",15,[[]]],[11,"drop","","",15,[[]]],[11,"vzip","","",15,[[]]],[11,"from","br::correct::exist","",3,[[]]],[11,"into","","",3,[[]]],[11,"borrow","","",3,[[]]],[11,"borrow_mut","","",3,[[]]],[11,"try_from","","",3,[[],["result",4]]],[11,"try_into","","",3,[[],["result",4]]],[11,"type_id","","",3,[[],["typeid",3]]],[11,"init","","",3,[[]]],[11,"deref","","",3,[[]]],[11,"deref_mut","","",3,[[]]],[11,"drop","","",3,[[]]],[11,"vzip","","",3,[[]]],[11,"from","br::correct::gap_size","",4,[[]]],[11,"into","","",4,[[]]],[11,"borrow","","",4,[[]]],[11,"borrow_mut","","",4,[[]]],[11,"try_from","","",4,[[],["result",4]]],[11,"try_into","","",4,[[],["result",4]]],[11,"type_id","","",4,[[],["typeid",3]]],[11,"init","","",4,[[]]],[11,"deref","","",4,[[]]],[11,"deref_mut","","",4,[[]]],[11,"drop","","",4,[[]]],[11,"vzip","","",4,[[]]],[11,"from","br::correct::graph","",5,[[]]],[11,"into","","",5,[[]]],[11,"borrow","","",5,[[]]],[11,"borrow_mut","","",5,[[]]],[11,"try_from","","",5,[[],["result",4]]],[11,"try_into","","",5,[[],["result",4]]],[11,"type_id","","",5,[[],["typeid",3]]],[11,"init","","",5,[[]]],[11,"deref","","",5,[[]]],[11,"deref_mut","","",5,[[]]],[11,"drop","","",5,[[]]],[11,"vzip","","",5,[[]]],[11,"from","br::correct::greedy","",6,[[]]],[11,"into","","",6,[[]]],[11,"borrow","","",6,[[]]],[11,"borrow_mut","","",6,[[]]],[11,"try_from","","",6,[[],["result",4]]],[11,"try_into","","",6,[[],["result",4]]],[11,"type_id","","",6,[[],["typeid",3]]],[11,"init","","",6,[[]]],[11,"deref","","",6,[[]]],[11,"deref_mut","","",6,[[]]],[11,"drop","","",6,[[]]],[11,"vzip","","",6,[[]]],[11,"from","br::error","",8,[[]]],[11,"into","","",8,[[]]],[11,"to_string","","",8,[[],["string",3]]],[11,"borrow","","",8,[[]]],[11,"borrow_mut","","",8,[[]]],[11,"try_from","","",8,[[],["result",4]]],[11,"try_into","","",8,[[],["result",4]]],[11,"type_id","","",8,[[],["typeid",3]]],[11,"init","","",8,[[]]],[11,"deref","","",8,[[]]],[11,"deref_mut","","",8,[[]]],[11,"drop","","",8,[[]]],[11,"as_error_source","","",8,[[],["error",8]]],[11,"vzip","","",8,[[]]],[11,"from","","",9,[[]]],[11,"into","","",9,[[]]],[11,"to_string","","",9,[[],["string",3]]],[11,"borrow","","",9,[[]]],[11,"borrow_mut","","",9,[[]]],[11,"try_from","","",9,[[],["result",4]]],[11,"try_into","","",9,[[],["result",4]]],[11,"type_id","","",9,[[],["typeid",3]]],[11,"init","","",9,[[]]],[11,"deref","","",9,[[]]],[11,"deref_mut","","",9,[[]]],[11,"drop","","",9,[[]]],[11,"as_error_source","","",9,[[],["error",8]]],[11,"vzip","","",9,[[]]],[11,"from","","",10,[[]]],[11,"into","","",10,[[]]],[11,"to_string","","",10,[[],["string",3]]],[11,"borrow","","",10,[[]]],[11,"borrow_mut","","",10,[[]]],[11,"try_from","","",10,[[],["result",4]]],[11,"try_into","","",10,[[],["result",4]]],[11,"type_id","","",10,[[],["typeid",3]]],[11,"init","","",10,[[]]],[11,"deref","","",10,[[]]],[11,"deref_mut","","",10,[[]]],[11,"drop","","",10,[[]]],[11,"as_error_source","","",10,[[],["error",8]]],[11,"vzip","","",10,[[]]],[11,"from","br::set::hash","",11,[[]]],[11,"into","","",11,[[]]],[11,"borrow","","",11,[[]]],[11,"borrow_mut","","",11,[[]]],[11,"try_from","","",11,[[],["result",4]]],[11,"try_into","","",11,[[],["result",4]]],[11,"type_id","","",11,[[],["typeid",3]]],[11,"init","","",11,[[]]],[11,"deref","","",11,[[]]],[11,"deref_mut","","",11,[[]]],[11,"drop","","",11,[[]]],[11,"vzip","","",11,[[]]],[11,"from","br::set::pcon","",12,[[]]],[11,"into","","",12,[[]]],[11,"borrow","","",12,[[]]],[11,"borrow_mut","","",12,[[]]],[11,"try_from","","",12,[[],["result",4]]],[11,"try_into","","",12,[[],["result",4]]],[11,"type_id","","",12,[[],["typeid",3]]],[11,"init","","",12,[[]]],[11,"deref","","",12,[[]]],[11,"deref_mut","","",12,[[]]],[11,"drop","","",12,[[]]],[11,"vzip","","",12,[[]]],[11,"valid_kmer","br::correct::exist","",3,[[],["boxkmerset",6]]],[11,"correct_error","","",3,[[],["option",4]]],[11,"valid_kmer","br::correct::gap_size","",4,[[],["boxkmerset",6]]],[11,"correct_error","","",4,[[],["option",4]]],[11,"valid_kmer","br::correct::graph","",5,[[],["boxkmerset",6]]],[11,"correct_error","","",5,[[],["option",4]]],[11,"k","br::correct::greedy","",6,[[]]],[11,"valid_kmer","","",6,[[],["boxkmerset",6]]],[11,"correct_error","","",6,[[],["option",4]]],[11,"init","br::correct::exist::one","",0,[[]]],[11,"c","","",0,[[]]],[11,"apply","","",0,[[["boxkmerset",6]],["option",4]]],[11,"correct","","",0,[[["boxkmerset",6]]]],[11,"init","br::correct::exist::two","",1,[[]]],[11,"c","","",1,[[]]],[11,"apply","","",1,[[["boxkmerset",6]],["option",4]]],[11,"correct","","",1,[[["boxkmerset",6]]]],[11,"get","br::set::hash","",11,[[]]],[11,"k","","",11,[[]]],[11,"get","br::set::pcon","",12,[[]]],[11,"k","","",12,[[]]],[11,"from","br::error","",8,[[["cli",4]]]],[11,"from","","",8,[[["io",4]]]],[11,"next_back","br::correct::exist::one","",14,[[],["option",4]]],[11,"next_back","br::correct::exist::two","",15,[[],["option",4]]],[11,"len","br::correct::exist::one","",14,[[]]],[11,"len","br::correct::exist::two","",15,[[]]],[11,"next","br::correct::exist::one","",14,[[],["option",4]]],[11,"size_hint","","",14,[[]]],[11,"nth","","",14,[[],["option",4]]],[11,"next","br::correct::exist::two","",15,[[],["option",4]]],[11,"size_hint","","",15,[[]]],[11,"nth","","",15,[[],["option",4]]],[11,"clone","br::correct::exist::one","",14,[[],["scenariooneiter",3]]],[11,"clone","","",0,[[],["scenarioone",4]]],[11,"clone","br::correct::exist::two","",15,[[],["scenariotwoiter",3]]],[11,"clone","","",1,[[],["scenariotwo",4]]],[11,"fmt","br::correct::exist::one","",0,[[["formatter",3]],["result",6]]],[11,"fmt","br::correct::exist::two","",1,[[["formatter",3]],["result",6]]],[11,"fmt","br::error","",8,[[["formatter",3]],["result",6]]],[11,"fmt","","",9,[[["formatter",3]],["result",6]]],[11,"fmt","","",10,[[["formatter",3]],["result",6]]],[11,"fmt","","",8,[[["formatter",3]],["result",6]]],[11,"fmt","","",9,[[["formatter",3]],["result",6]]],[11,"fmt","","",10,[[["formatter",3]],["result",6]]],[11,"source","","",8,[[],[["error",8],["option",4]]]],[11,"iter","br::correct::exist::one","",0,[[],["scenariooneiter",3]]],[11,"iter","br::correct::exist::two","",1,[[],["scenariotwoiter",3]]]],"p":[[4,"ScenarioOne"],[4,"ScenarioTwo"],[8,"Scenario"],[3,"Exist"],[3,"GapSize"],[3,"Graph"],[3,"Greedy"],[8,"Corrector"],[4,"Error"],[4,"Cli"],[4,"IO"],[3,"Hash"],[3,"Pcon"],[8,"KmerSet"],[3,"ScenarioOneIter"],[3,"ScenarioTwoIter"]]},\
"br_large":{"doc":"","i":[[3,"Command","br_large","",null,null],[12,"inputs","","",0,null],[12,"outputs","","",0,null],[12,"kmer_solid","","",0,null],[12,"kmer_size","","",0,null],[12,"methods","","",0,null],[12,"confirm","","",0,null],[12,"max_search","","",0,null],[12,"two_side","","",0,null],[12,"threads","","",0,null],[12,"record_buffer","","",0,null],[12,"verbosity","","",0,null],[5,"main","","",null,[[],["result",6]]],[11,"from","","",0,[[]]],[11,"into","","",0,[[]]],[11,"borrow","","",0,[[]]],[11,"borrow_mut","","",0,[[]]],[11,"try_from","","",0,[[],["result",4]]],[11,"try_into","","",0,[[],["result",4]]],[11,"type_id","","",0,[[],["typeid",3]]],[11,"init","","",0,[[]]],[11,"deref","","",0,[[]]],[11,"deref_mut","","",0,[[]]],[11,"drop","","",0,[[]]],[11,"vzip","","",0,[[]]],[11,"fmt","","",0,[[["formatter",3]],["result",6]]],[11,"into_app","","",0,[[],["app",3]]],[11,"augment_clap","","",0,[[["app",3]],["app",3]]],[11,"from_arg_matches","","",0,[[["argmatches",3]]]]],"p":[[3,"Command"]]}\
}');
addSearchOptions(searchIndex);initSearch(searchIndex);