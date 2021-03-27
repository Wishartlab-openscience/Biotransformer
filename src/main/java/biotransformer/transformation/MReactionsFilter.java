/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.transformation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import org._3pq.jgrapht.alg.ConnectivityInspector;
import org._3pq.jgrapht.graph.SimpleDirectedGraph;

import biotransformer.biosystems.BioSystem;
import biotransformer.biosystems.BioSystem.BioSystemName;
//import biotransformer.transformation.MRPatterns.ReactionName;


public class MReactionsFilter {

	public LinkedHashMap<String, Object> btRelativePrecedenceRules;
	public ArrayList<String> btStrictPrecedences;
	public BioSystemName bioSysName;
	
	public MReactionsFilter(BioSystem bioSys){
		this.bioSysName = bioSys.name;
		this.setPriorityHash(bioSys);
	}
	
	
	@SuppressWarnings("unchecked")
	public void setPriorityHash(BioSystem bioSys){
			
		this.btRelativePrecedenceRules = (LinkedHashMap<String, Object>) bioSys.reactionPrecedenceRules.get("relative");
		this.btStrictPrecedences = (ArrayList<String>) bioSys.reactionPrecedenceRules.get("strict");
		
//		System.out.println(this.btRelativePrecedenceRules);
//		System.out.println(this.btRelativePrecedenceRules.size());
	}
	
	
	private static LinkedHashMap<String, String[]> unfoldReactionPriorityMap(LinkedHashMap<String, Object> rpmap){
		LinkedHashMap<String, String[]> unfolded = new LinkedHashMap<String, String[]>();
		
		// Following the example on this page: https://github.com/jgrapht/jgrapht/wiki/DirectedGraphDemo
		SimpleDirectedGraph sdg = new SimpleDirectedGraph();
		
		for(String rn : rpmap.keySet()){
			String[] rn_children = (String[]) ((LinkedHashMap<String, Object>) rpmap.get(rn)).get("reactions");
			if(!sdg.containsVertex(rn)){
				sdg.addVertex(rn);
			}
			
			for(int i = 0; i < rn_children.length; i++){
				if(!sdg.containsVertex(rn_children[i])){
					sdg.addVertex(rn_children[i]);
				}
				sdg.addEdge(rn, rn_children[i]);
			}			
		}		
		System.out.println(sdg.edgeSet().size() + " edges from " + sdg.vertexSet().size() + " vertices.");
	
		ConnectivityInspector ci = new ConnectivityInspector(sdg);
	
		for( String rn : rpmap.keySet() ){
			Set<String> des =  new  HashSet<String>();
			for(String rn2: rpmap.keySet()){
				if(rn != rn2 && ci.pathExists(rn, rn2)){
					des.add(rn2);
				}
				unfolded.put(rn, des.toArray(new String[des.size()]));				
			}
//			System.out.println(rn + " = " + unfolded.get(rn).length);			
		}
//		System.out.println(rpmap.size());

		return unfolded;
	}
	

	
	public LinkedHashMap<String, MetabolicReaction> filterReactions(ArrayList<MetabolicReaction> reactions){
		return filterReactions(reactions, this.btRelativePrecedenceRules, this.btStrictPrecedences);
	}

	public static LinkedHashMap<String, MetabolicReaction> filterReactions(ArrayList<MetabolicReaction> reactions, LinkedHashMap<String, Object> relativeRules, ArrayList<String> strictPrecedences ){

		LinkedHashMap<String, MetabolicReaction> reactionsHash = new LinkedHashMap<String, MetabolicReaction>();
		LinkedHashMap<String, MetabolicReaction> filteredReactions = new LinkedHashMap<String, MetabolicReaction>();
		ArrayList<String> rNames = new ArrayList<String>();
		
		for(MetabolicReaction i : reactions){
			reactionsHash.put(String.valueOf(i.name), i);
			rNames.add(String.valueOf(i.name));
		}

		ArrayList<String> rNamesWithUMStrictReasoning = new ArrayList<String>();
		
		for(String m : strictPrecedences) {
			if(rNames.contains(m)) {
				rNamesWithUMStrictReasoning.add(m);
			}
		}

		
//		System.err.println("rNames : " + rNames);
		
		if(rNamesWithUMStrictReasoning.size()>0){
			LinkedHashMap<String, MetabolicReaction> filteredReactionsStrict = new LinkedHashMap<String, MetabolicReaction>();
			for(String rn : rNamesWithUMStrictReasoning){
				filteredReactionsStrict.put(rn, reactionsHash.get(rn));
			}
//			System.out.println("rNamesWithUMStrictReasoning : " + rNamesWithUMStrictReasoning);
			return filteredReactionsStrict;
		} else {
			
			ArrayList<String> parentNodes = new ArrayList<String>();
			ArrayList<String> childNodes = new ArrayList<String>();
			ArrayList<String> nodesWithNoParentOrChild = new ArrayList<String>();

			for(String rN : rNames){
				
				if(relativeRules.containsKey(rN)){
//				System.out.println(rN);
//				System.out.println("REACTIONS=" + ((LinkedHashMap<String, Object>) relativeRules.get(rN)).get("reactions"));
				ArrayList<String> rNpriorities = (ArrayList<String>) ((LinkedHashMap<String, Object>) relativeRules.get(rN)).get("reactions");
				
				
					for(String rN2 : rNames){
						if(rN2 != rN){
							if(rNpriorities.contains(rN2)){
								parentNodes.add(rN);
							} else {
								
								if(relativeRules.containsKey(rN2)) {
									ArrayList<String> rN2priorities = (ArrayList<String>) ((LinkedHashMap<String, Object>) relativeRules.get(rN2)).get("reactions");
									if (rN2priorities.contains(rN))	{
										childNodes.add(rN);
									}									
								}
								
							}
						}
					}
				
				}
				if(!(parentNodes.contains(rN) || childNodes.contains(rN))){
					nodesWithNoParentOrChild.add(rN);
				}
			}
			
//			System.out.println("parentNodes : " + parentNodes);
//			System.out.println("childNodes : " + childNodes);
//			System.out.println("nodesWithNoParentOrChild : " + nodesWithNoParentOrChild);
			
			ArrayList<String> rNames_reduced  = (ArrayList<String>) rNames.clone();
			
			// This will avoid remove nodes Xi that are children of Yi, but not children of any other reaction.
			for(String r: rNames){
				if(parentNodes.contains(r) && childNodes.contains(r)){
					rNames_reduced.remove(r);
				}
			}
			
			
			ArrayList<String> toRemove = new ArrayList<String>();
			for(String r : rNames){
				for(String r2 : rNames){
					if(r2 != r){
						if(relativeRules.containsKey(r)) {
							ArrayList<String> rpriorities = (ArrayList<String>) ((LinkedHashMap<String, Object>) relativeRules.get(r)).get("reactions");						
							if(rpriorities.contains(r2)){
								toRemove.add(r2);
							}
						}

					}
				}
			}
				
//			System.err.println("To remove");
//			for(ReactionName rr : toRemove){
//				System.err.println("\t" + rr);
//			}
			
			rNames_reduced.removeAll(toRemove);
			
//			System.out.println("toRemove : " + toRemove);
//			System.out.println("rNames_reduced : " + rNames_reduced);
			
			for( String rs : rNames_reduced){
				filteredReactions.put(rs, reactionsHash.get(rs));
			}
			
		}		
//		System.out.println("filteredReactions : " + filteredReactions);
		return filteredReactions;
	}

	
	
//	public static LinkedHashMap<ReactionName, MetabolicReaction> filterReactions(ArrayList<MetabolicReaction> reactions, LinkedHashMap<ReactionName, ReactionName[]> btph ){
//
//		LinkedHashMap<ReactionName, MetabolicReaction> filteredReactions = new LinkedHashMap<ReactionName, MetabolicReaction>();
//		ArrayList<ReactionName> rNames = new ArrayList<ReactionName>();
//		
//		for(MetabolicReaction i : reactions){
//			filteredReactions.put(ReactionName.valueOf(i.name), i);
//			rNames.add(ReactionName.valueOf(i.name));
//		}
//		
//		
//		ArrayList<ReactionName> rNamesWithUMStrictReasoning = new ArrayList<ReactionName>();
//		if(rNames.contains(ReactionName.EAWAG_RULE_BT0003)){
//			rNamesWithUMStrictReasoning.add(ReactionName.EAWAG_RULE_BT0003);
//		}
//		if(rNames.contains(ReactionName.EAWAG_RULE_BT0255_PATTERN1)){
//			rNamesWithUMStrictReasoning.add(ReactionName.EAWAG_RULE_BT0255_PATTERN1);
//		}
//		if(rNames.contains(ReactionName.EAWAG_RULE_BT0255_PATTERN2)){
//			rNamesWithUMStrictReasoning.add(ReactionName.EAWAG_RULE_BT0255_PATTERN2);
//		}		
//		if(rNamesWithUMStrictReasoning.size()>0){
//			LinkedHashMap<ReactionName, MetabolicReaction> filteredReactionsStrict = new LinkedHashMap<ReactionName, MetabolicReaction>();
//			for(ReactionName rn : rNamesWithUMStrictReasoning){
//				filteredReactionsStrict.put(rn, filteredReactions.get(rn));
//			}
//			
//			return filteredReactionsStrict;
//		} else {
//		
//			
//			
//			
//			ArrayList<ReactionName> toRemove = new ArrayList<ReactionName>();
//			for(ReactionName rN : rNames){
//				System.out.println(rN);
//				if(btph.get(rN) !=null ){
//					for(int n = 0; n < btph.get(rN).length; n++){
//						System.out.println("\t" + btph.get(rN)[n]);
//						toRemove.add( btph.get(rN)[n] );
//					}					
//				}
//
//			}
//			
//			for(ReactionName r : toRemove){
//				filteredReactions.remove(r);
//			}
//			
//		}
//		
//		
//		return filteredReactions;
//	}
	
	
	
//	public static LinkedHashMap<ReactionName, MetabolicReaction> filterReactions(ArrayList<MetabolicReaction> reactions, LinkedHashMap<ReactionName, ReactionName[]> btph ){
//
//		LinkedHashMap<ReactionName, MetabolicReaction> filteredReactions = new LinkedHashMap<ReactionName, MetabolicReaction>();
//		ArrayList<ReactionName> rNames = new ArrayList<ReactionName>();
//		
//		for(MetabolicReaction i : reactions){
//			filteredReactions.put(ReactionName.valueOf(i.name), i);
//			rNames.add(ReactionName.valueOf(i.name));
//		}
//
//		ArrayList<ReactionName> rNamesWithUMStrictReasoning = new ArrayList<ReactionName>();
//		if(rNames.contains(ReactionName.EAWAG_RULE_BT0003)){
//			rNamesWithUMStrictReasoning.add(ReactionName.EAWAG_RULE_BT0003);
//		}
//		if(rNames.contains(ReactionName.EAWAG_RULE_BT0255_PATTERN1)){
//			rNamesWithUMStrictReasoning.add(ReactionName.EAWAG_RULE_BT0255_PATTERN1);
//		}
//		if(rNames.contains(ReactionName.EAWAG_RULE_BT0255_PATTERN2)){
//			rNamesWithUMStrictReasoning.add(ReactionName.EAWAG_RULE_BT0255_PATTERN2);
//		}		
//		if(rNamesWithUMStrictReasoning.size()>0){
//			LinkedHashMap<ReactionName, MetabolicReaction> filteredReactionsStrict = new LinkedHashMap<ReactionName, MetabolicReaction>();
//			for(ReactionName rn : rNamesWithUMStrictReasoning){
//				filteredReactionsStrict.put(rn, filteredReactions.get(rn));
//			}
//			
//			return filteredReactionsStrict;
//		} else {
//		
//			int count = reactions.size();
//			int index = 0;
////			System.err.println("FILTERING....\n");
//				while(index< reactions.size() && !rNames.isEmpty()){
////					System.out.println("Reaction: " + rNames.get(index));
//					if(btph.containsKey( rNames.get(index) )){
//						List<ReactionName> current = Arrays.asList(btph.get( rNames.get(index)));			
////						System.out.println("# of reactions with lower priority: " + current.size());
//						for( int k = 0; k < current.size(); k++){
//							
//							if(btph.containsKey(current.get(k))){
////								System.out.println(current.get(k) + " dominates over " + btph.get(current.get(k)).length + " reactions");
//								for(ReactionName r : btph.get(current.get(k))){
//									if(filteredReactions.containsKey(r)){
////										System.err.println("REMOVING REACTION: " + r);
//										filteredReactions.remove(r);
//										count--;										
//									}
//									
//
//								}
//							}
//							filteredReactions.remove(current.get(k));
//							count--;
//						}
//					}
//					index++;
//				}
////				System.err.println("\nRetaining " + filteredReactions.size() + " reactions, with filter set to true");
////				for(ReactionName nn : filteredReactions.keySet()){
////					System.out.println("Retained " + nn);
////				}
//				
//				return filteredReactions;
//		}
//		
//	}
//	
	public LinkedHashMap<String, MetabolicReaction> filterReactions(LinkedHashMap<String, MetabolicReaction> reactions){
		return filterReactions(reactions, unfoldReactionPriorityMap(this.btRelativePrecedenceRules));
	}
	
	public static LinkedHashMap<String, MetabolicReaction> filterReactions(LinkedHashMap<String, MetabolicReaction> reactions, LinkedHashMap<String, String[]> btph){
		LinkedHashMap<String, MetabolicReaction> filteredReactions = new LinkedHashMap<String, MetabolicReaction>();

		ArrayList<String> rNames = (ArrayList<String>) reactions.keySet();
		List<String> rNamesFiltered = new ArrayList<String>(reactions.keySet());

		int count = reactions.size();
		int index = 0;
		
		while(index< reactions.size() && !rNames.isEmpty()){	
			if(btph.containsKey( rNames.get(index))){
			List<String> current = Arrays.asList(btph.get( rNames.get(index)));
//			System.out.println("# of reactions with lower priority: " + current.size());	
			for( int k = 0; k < current.size(); k++){
//					System.err.println("REMOVING REACTION: " + current.get(k));
				
				if(btph.containsKey(current.get(k))){
					for(String r : btph.get(current.get(k))){
						filteredReactions.remove(r);
						count--;
					}
				}
				
				rNamesFiltered.remove(current.get(k));
					
					count--;
				}
			}
			index++;
		}

		for (int z = 0; z < rNamesFiltered.size(); z++){
			filteredReactions.put(rNames.get(z), reactions.get(rNames.get(z)));
		}
//		System.err.println("Retaining " + filteredReactions.size() + " reactions, with filter set to true");
		return filteredReactions;
		
	}

	public LinkedHashMap<String, MetabolicReaction> filterReactionsLSM(LinkedHashMap<String, MetabolicReaction> reactions){
		LinkedHashMap<String, MetabolicReaction> filteredReactions = new LinkedHashMap<String, MetabolicReaction>();
		
		List<String> rNames = new ArrayList<String>(reactions.keySet());
		List<String> rNamesFiltered = new ArrayList<String>(reactions.keySet());

		int count = reactions.size();
		int index = 0;
		
		while(index< reactions.size() && !rNames.isEmpty()){
			System.err.println("rNamesFiltered: " + rNamesFiltered.size());
			System.err.println("Index: " + rNames.get(index));
			if(this.btRelativePrecedenceRules.containsKey(String.valueOf(rNames.get(index)))){
//				System.err.println(rNames.get(index) + " dominates over reactions in btPriorityHash");
				List<String> current = 
					Arrays.asList((String[])this.btRelativePrecedenceRules.get(String.valueOf(rNames.get(index))));
				for( int k = 0; k < current.size(); k++){
//					System.err.println("REMOVING REACTION: " + current.get(k));
					rNamesFiltered.remove(current.get(k).toString());
					count--;
				}		
			}		
			index++;
		}

		for (int z = 0; z < rNamesFiltered.size(); z++){
			filteredReactions.put(String.valueOf(rNames.get(z)), reactions.get(rNames.get(z)));
		}
//		System.err.println("Applied " + filteredReactions.size() + " reactions, with filter set to true");
		return filteredReactions;	
	}

	
	public MReactionsFilter() {
		// TODO Auto-generated constructor stub
	}

}
