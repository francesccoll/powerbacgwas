#!/usr/bin/perl -w

# This Perl script makes use of Bio::TreeIO library to parse, annotate internal nodes (replacing bootstrapping values) and write the resulting tree. A file describing nodes will be produced to. Note that this script has been tested on trees of newick format produced by fasttree.

use strict;
use Bio::TreeIO;


my $argn = $#ARGV+1;

if($argn != 3)
{
	print "\nperl annotate_nodes_newick.pl input_tree node_prefix output_node_file\n";
	print "\nDESCRIPTION:";
	print "\nThis Perl script makes use of Bio::TreeIO library to parse, annotate internal nodes (replacing bootstrapping values) and write the resulting tree. A file describing nodes will be produced to. Note that this script has been tested on trees of newick format produced by fasttree.";
	print "\n\nARGUMENTS:";
	print "\n\tinput_tree\t\tInput tree of newick format";
	print "\n\tnode_prefix\t\tNode prefix used to label internal nodes (e.g. \"node\")";
	print "\n\toutput_node_file\tOutput node file";
	print "\n\n";
	exit(0);
}
my $input_tree = $ARGV[0];
my $node_prefix = $ARGV[1];
my $output_node_file = $ARGV[2];
chomp($output_node_file);
unless(-e $input_tree){ print $input_tree." could not be found\n"; exit(0); }

my %lhBootstrap = (); 
 
# Reading Input Tree
my $treeio = new Bio::TreeIO(-file   => $input_tree,
                            -format => "newick");
my $tree = $treeio->next_tree;

my @nodes  = $tree->get_nodes; 		# array of Bio::Tree::NodeI objects
print "Number of nodes: ".($#nodes+1)."\n";
my @leaves = $tree->get_leaf_nodes;	# Returns the leaves (tips) of the tree, Array of Bio::Tree::NodeI objects
print "Number of leaf nodes: ".($#leaves+1)."\n";
my $root   = $tree->get_root_node;
my $root_id = "Nonde";
if(defined $root->id){ $root_id = $root->id; }
print "Root node ID: ".$root_id."\n";
my $size = $tree->number_nodes;		# Find the number of nodes in the tree (int)
print "Number of nodes: ".$size."\n";
my $size2 = $tree->total_branch_length;	# Returns the sum of the length of all branches (int)
print "Sum of the length of all branches ".$size2."\n";
my $tree_id = "None";
if(defined $tree->id()){ $tree_id = $tree->id(); }
print "Id value for the tree ".$tree_id."\n";





$treeio = new Bio::TreeIO(-file   => $input_tree,
                            -format => "newick");

my $out = new Bio::TreeIO(-file => '>outtre',
                          -format => 'newick');

my $count = 0;
while($tree = $treeio->next_tree) 
{
 my @nodes = $tree->get_nodes;
 for my $node ( $tree->get_nodes )
 {
	my $node_id	= $node->id;
	my $boostrap	= 0;
	if(defined $node_id) { $boostrap = $node_id; }
	
	unless($node->is_Leaf)
	{ 
		#if(defined $node_id)
		#{
			$count++;
			my $new_id = $node_prefix.$count;
			print "Node ".$count." ID:".$node_id." --> ";
			$node->id($new_id);
			$node->bootstrap($boostrap);
			print $node->id." --> ".$node->bootstrap."\n";
			$lhBootstrap{$new_id} = $boostrap;
		#}
	}
 }
	$out->write_tree($tree);
}


# Saving new tree and internal node attributes


$treeio = new Bio::TreeIO(-file   => 'outtre',
                            -format => "newick");

$output_node_file = $output_node_file.".txt";

open O, ">$output_node_file";
my $header = "node_Id\tBoostrap_value\tDepth\tNum_ancestor\tNum_descendents\tNum_descendents_leaves\tAncetral_nodes_Id\tDescendent_nodes_Id\tDescendent_nodes_leaves_Id\n";
print O $header;


while($tree = $treeio->next_tree) 
{
 for my $node ( $tree->get_nodes )
 {
	my $branch_len 	= $node->branch_length;
	my $node_id	= $node->id;

	unless($node->is_Leaf)
	{ 
			my @nodes3 = $node->get_all_Descendents; # Recursively fetch all the nodes and their descendents
			# Gather descendent node IDs
			my $des_nodes_ids = "";
			my $des_nodes_leaves_ids = "";
			my $des_nodes_leaves_num = 0;
			foreach my $des_node (@nodes3)
			{
				$des_nodes_ids = $des_nodes_ids.$des_node->id.";";
				if($des_node->is_Leaf)
				{
					$des_nodes_leaves_ids = $des_nodes_leaves_ids.$des_node->id.";";
					$des_nodes_leaves_num++;
				}
			}
			$des_nodes_ids =~ s/\;$//g; # removing last ;
			$des_nodes_leaves_ids =~ s/\;$//g;
		
			my $anc_node = $node->ancestor;
			my $anc_node_id = "None";
			if(defined $anc_node){ $anc_node_id = $anc_node->id; }
			my $des_count = $node->descendent_count; # Counts the number of descendents a node has
			my $len = $node->depth; # Depth is the distance from this node to the root.
			my $node_bst = "0";
			if(defined $node->bootstrap){ $node_bst = $node->bootstrap; }
			print "\$node_id: ".$node_id."\n";

			my $line_node_inf = $node_id."\t".$lhBootstrap{$node_id}."\t".$len."\t"."1"."\t".$des_count."\t".$des_nodes_leaves_num."\t".$anc_node_id."\t".$des_nodes_ids."\t".$des_nodes_leaves_ids."\n";
			print $line_node_inf;
			print O $line_node_inf;
	}
 }
}

close O;


# Node names do not match SNPSTI node values (assigned sequentially as they come)








# http://doc.bioperl.org/releases/bioperl-1.4/Bio/Tree/TreeFunctionsI.html

# find all the nodes named 'node1' (there should be only one)
#my @nodes = $tree->find_node(-id => 'node1');
# find all the nodes which have description 'BMP'
#my @nodes = $tree->find_node(-description => 'BMP');
# find all the nodes with bootstrap value of 70
#my @nodes = $tree->find_node(-bootstrap => 70);








