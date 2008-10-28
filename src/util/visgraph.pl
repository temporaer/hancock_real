use GraphViz;
use List::Util qw/reduce/;
use Data::Dumper;

$adj  = `head -n 1 ../../build/cb/out.dat`;
$path = `head -n 2 ../../build/cb/out.dat | tail -n 1`;
$path =~ s/^.*\[\s*//;
$path =~ s/[\];]//g;
print "$path\n";

sub getmat{
	$_ = shift();
	s/(\s*;\s*)+$//g;
	s/^.*\[//g;
	s/^/[/;
	s/(\d)\s+(\d)/\1,\2/g;
	s/(\d)\s+(\d)/\1,\2/g;
	s/;/],[/g;
	return $_;
}

#create your adjacency matrix, nodes do whatever else you need.
my @admat = (
	eval(getmat($adj))
);

@nodes = map{ $_ } (0..$#admat);

#create a new graph
my $graph = GraphViz->new( 
	layout => 'circo',
	directed => 0);

#add all your nodes to the graph
foreach (@nodes) {
    $graph->add_node($_);
}

#create edges as defined by your adjacecncy matrix
for($i=0; $i<@nodes; $i++){
    for($j=$i+1; $j<@nodes; $j++){
        if($admat[$i][$j] eq 1){
            print "adding edge from $i to $j\n";
			$graph->add_edge($nodes[$i] => $nodes[$j], color => 'gray');
		}
    }
}


print "Create paths...\n";
drawpath('red',  split(/\s+/, $path));

sub drawpath{
	$color = shift;
	print "  drawing $color path...\n";
	reduce{ 
		if($admat[$a][$b] eq 1){
			$graph->add_edge($a => $b, color => $color, arrowhead => 'normal'); 
		}else{
			print "    jumping from ".($a)." to ".($b)."\n";
		};
		$b
	} @_;
};


#render the graph to a png-file
$graph->as_png("pretty.png");
