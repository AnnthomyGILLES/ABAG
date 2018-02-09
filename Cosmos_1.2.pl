########################################################################
#			            		COSMOS 1.2		                       #		
########################################################################

use strict;
use warnings;no warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Transcript;
use List::Util 'first'; 
use LWP::UserAgent;
use LWP::Simple qw( $ua );
#use LWP::Debug qw(+); #Le chargement de ce module force LWP::UserAgent à afficher ce qu'il fait. 

use List::Uniq ':all';

use Bio::DB::EUtilities;

#-----------------Ouverture en lecture du fichier----------------------#
open (file_input,'<Genesymbols.txt')or die "Probleme pour ouvrir: $!";
open(FILE_RESULT, '>results.html');


#---------------Creation de l'entête du fichier HTML-------------------#

my $header='<html>
<head>
<meta http-equiv="Content-type" content="text/html; charset=utf-8">
		<title>Cosmos 1.2</title>

 
</head>

<body>
	<h1  style="text-align: center;"> Cosmos 1.2</h1>
<div style="text-align: center;height: 100px"><a href="http://masterbioinfo.formations.univ-rouen.fr/">
<img src="rsz_2logo.png" class="imageflottante" alt="Logo Master BIM" title="Visiter le site du Master BIM"></a></div>
<table id="example" class="display" cellspacing="0" width="100%">
<h4  style="text-align: center;"> Cosmos, Less Is More...</h1>

<thead>
	<tr>';
#---------------Creation de l'entête du tableau de sortie (en HTML)-------------------#

my @headers=("Gene Symbol (NCBI)","Organism","Gene ID (NCBI)","Official full name","ENSG_","Genome Browser","ENST_","ENSP_","refseq NM_","Refseq NP_", "Uniprot ID", "Protein Name","PDB ID","STRING","Gene ontology(cellular component)", "Gene ontology(cellular Function)", "Gene ontology (biological process)");

print FILE_RESULT "$header\n";


	foreach my $header_name(@headers){#Boucle qui parcourt chaque élément du tableau d'entête
		print FILE_RESULT "<th>$header_name</th>\n";#Ecriture des éléments du tableau entre des balises <th> afin de créer l'entête du tableau
	}

print FILE_RESULT "
</tr>
	</thead>
	<tbody>	
";

#---------------Déclaration des variables et tableaux-------------------#

######Variables Ensembl#####
my ($ensg,$enst,$ensp,$transcript,$translation,$gen,$trans,$gene_symbol,$genebrows_url);
my (@xrefrna,@xrefprot,$xrefrna,$xrefprot);
my ($response,$gene,$gstring);
my $file;
######Variables PDB#####
my ($list_pdbID);
########Variables NCBI#####
my ($refseq_mrna,$refseq_prot,$gene_id,$id,$off_full_name);
########Variables Uniprot#####
my ($uniprot_source,$uniprot_ids,$protein_name_unip,$id_unip,@protein_name,$uniprot_acc,@uniprot_ID);
########Variables STRING#####
my (@tab_string_url,$string_url);
########Variables QuickGO#####
my ($tab_quickgo_process,$tab_quickgo_function,$tab_quickgo_component,$tab_quickgo_aspect);


#---------------Parcours du fichier-------------------#
while($file=<file_input>){
	print "Recuperation du gene symbol et de l'organisme dans le fichier... \n";
	next if( $file =~ /^(\s)*$/);#Si il y a une ligne vide, fin de la boucle.
	my @cols = split(/\t/, $file); #eclate la chaîne de caractère en sous parties, et la place dans un tableau 
	my $species = $cols[1];chomp $species; #Recupère le second element du tableau (organisme) et le place dans la variable correspondante
	
	my $source  = 'core'; # core or vega
	my $gene_symbol=$cols[0];chomp $gene_symbol;#Recupère le 1er element du tableau (gene symbole) et le place dans la variable correspondante
	print "le Gene symbol est: ".$gene_symbol." et l'organisme: ".$species."\n";

	
########################################################################
#			            Ensembl				                           #		
########################################################################

	my $r="Bio::EnsEMBL::Registry";# Cree un objet de type "Bio::EnsEMBL::Registry".
	#my $s="Bio::SeqIO";
	#Recupère le registre depuis la base de donnees
	$r->load_registry_from_db(-host=>"ensembldb.ensembl.org", #Nom de l'hôte de la base de donnees Ensembl
											-user=>"anonymous", #Nom d'utilisateur 
											-verbose=>"0");#0 afin de masquer les informations relatives à la connexions
	#--------------------------------------------------------------------------------------------#
	#Pour chaque objet, il existe un adaptateur.
	#L'adaptateur est lui-même un objet "intermediaire" qui contient,
	#entre autres, toutes les methodes capables de retourner un objet du type de l'adaptateur.
	#--------------------------------------------------------------------------------------------#
								
	#Recupère un adaptateur correspondant à l'objet qui nous interesse
	my $gene_adaptor   = $r->get_adaptor($species, $source,"Gene");#ici, un gène
	my $transcript_adaptor = $r->get_adaptor($species, $source,'Transcript' ); #ici, un transcrit.
	print "Extraction des IDentifiants EnsEMBL et reconstruction d'URL Genome Browser............ ";

	$gene=$gene_adaptor->fetch_by_display_label($gene_symbol);
	my @refseq_to_uniprot=();
	$ensg=$gene->stable_id;#Récupération de l'identifiant ENSG_
	my $seq_region = $gene->slice->seq_region_name();#Obtient la localisaiton chromosomique
	my $start = $gene->start();# Localisation: Region start
	my $end = $gene->end();#Localisation: Region end
	my $gstring="$seq_region:$start-$end"; #Obtient la localisation precise du gène sur le genome
	my $species_sub=$species;chomp $species_sub;
	$species_sub=~s/\s/_/; #Remplace les espaces par des _


	$genebrows_url="<a href=http://www.ensembl.org/$species_sub/Location/View?db=core;g=$ensg;r=$gstring>$gstring</a>";#Construit l'url vers Genome Browser

	my $ensg_canonical_transcript= $gene_adaptor->fetch_by_stable_id($ensg);
	my $canonical = $ensg_canonical_transcript->canonical_transcript();#Obtient le transcrit canonique
	my $canonical_transcript=$canonical->stable_id;	
	
	#Initialisation de 2 tableaux;
	@xrefrna=();@xrefprot=();
		foreach $trans (@{$gene->get_all_Transcripts}) { #Parcours de la liste de Transcripts
			$enst= $trans->stable_id;chomp $enst;#Extraction de l'ENST_
			$transcript = $transcript_adaptor->fetch_by_stable_id("$enst");#Prend en argument une variable et retourne un objet stable id de Bio::EnsEMBL::Transcript-->ENST_
			$translation = $transcript->translation;
			
	
			# Obtenir les ENST_ et ENSP_
			if(defined($translation)) { #Si la variable a une valeur
					$ensp=$translation->display_id();chomp $ensp;#Extraction de l'identifiant ENSP_
					#---------------Ajout par itération des ENST_ et ENSP_ dans un tableau-------------------#
					push(@xrefrna,"<a href=http://www.ensembl.org/Multi/Search/Results?q=$enst;site=ensembl;page=1;facet_feature_type=Transcript>$enst</a>"); #Reconstruit l'URL de l'ID ENST_ et ajoute la variable dans le tableau
					push(@refseq_to_uniprot,$ensp);
					push(@xrefprot,"<a href=http://www.ensembl.org/Multi/Search/Results?q=$ensp;site=ensembl>$ensp</a>");#Ajout de la variable (ENSP_) dans le tableau	
			}
			#else {
					#my $na="NA";
					#push(@xrefrna,$na);
					#push(@xrefprot,$na);
			#}		
		}
	print "Extraction terminee \n";

		
########################################################################
#							UNIPROT          						   #
########################################################################
	print "Extraction des IDentifiants Uniprot................";
	$xrefrna = join(' ', @xrefrna);#Joint les elements dans une chaîne avec un separateur ' '
	$xrefprot = join(' ', @xrefprot);#Idem
	my $ref_to_unip = join(' ', @refseq_to_uniprot);#Idem
		if($ref_to_unip !=""){#Si la variable n'est pas vide
			@uniprot_ID=();@tab_string_url=();$protein_name_unip=();
				foreach my $cle (@refseq_to_uniprot){
					chomp $cle;
					my $base = 'http://www.uniprot.org'; #La base qui sera interrogée
					my $tool = 'mapping'; #L'outils Mapping

					my $params = {
					from => 'ENSEMBL_PRO_ID',#De l'ID ENSP_
					to => 'ACC',#Au numero d'accession uniprot
					format => 'list',
					query => ''.$cle.''
					};

					my $contact = 'gotama@hotmail.com'; # Insertion de l'email
					my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");# initialisation de l'agent
					push @{$agent->requests_redirectable}, 'POST';

					my $response = $agent->post("$base/$tool/", $params);
						while (my $wait = $response->header('Retry-After')) {#Headers représentent les en-têtes de la réponse 
							print STDERR "Waiting ($wait)...\n";
							sleep $wait;
							$response = $agent->get($response->base);
						}

						unless( $response->is_success ) { #is_success permet de savoir si la requête s'est bien deroulee
							print "\n an error occurred: ", $response->status_line, "\n";#Dans le cas contraire, affichage d'un message d'erreur avec status_line
						}		
					$uniprot_acc= $response->content;
					push(@uniprot_ID,$uniprot_acc);
					
					
			}
		}
		else {#Si pour une raison quelconque l'Uniprot ID n'est pas recupere, le programme passe par une cross ref 
			my @xrefs=@{$gene->get_all_xrefs('Uniprot/S%')};#Les ID Uniprot sont places dans un tableau
			@tab_string_url=();@uniprot_ID=();#Vide le tableau avant d'entamer la bouche
				foreach my $cle (@xrefs){#$cle prend la valeur de chaque element de la liste @xrefs
					$uniprot_acc=$cle->display_id();chomp ($uniprot_acc);#l'ID Uniprot est extrait
					push(@uniprot_ID,$uniprot_acc);#Ajoute l'ID Uniprot dans un nouveau tableau
					@uniprot_ID=uniq(@uniprot_ID);#Supprime la redondance du tableau
					push(@tab_string_url,"<a href=http://string-db.org/api/image/network?identifier=$uniprot_acc>$uniprot_acc</a>");#Ajout de la variable (url) dans le tableau (Iteration)		
					@tab_string_url=uniq(@tab_string_url);

				}
		}
		@uniprot_ID=uniq(@uniprot_ID);#Suppression de la redondance du tableau
		$id_unip = join(' ', @uniprot_ID);#Joint les elements dtu tableau dans une chaîne avec un separateur ' '
		@uniprot_ID=split(/ \t/,$id_unip);
	print "Extraction terminee \n";

		
########################################################################
#							UNIPROT - Protein names	         		   #
########################################################################
	print "Extraction du nom de la proteine sur Uniprot........";
	#--------------Initialisation de l'agent------------------#
	my $ua = new LWP::UserAgent;
		
	#---------------Fabrication de la requête-----------------#
	my $response = $ua->get('http://www.uniprot.org/uniprot/?format=tab&query='.$canonical_transcript.'&columns=protein%20names');
	#---------------En cas d'echec de la requête--------------#

		unless ($response->is_success) {
			die $response->status_line;
		}
	#------------Attribution du contenu de la page dans $protein_name_unip-----#
	$protein_name_unip=$response->content;
	print "Extraction terminee \n";
	
########################################################################
#							     STRING 	         		           #
########################################################################
	print "Reconstruction d'URL STRING.................";

		foreach my $uniprot_acc (@uniprot_ID){#Parcours du tableau
			push(@tab_string_url,"<a href=http://string-db.org/api/image/network?identifier=$uniprot_acc>$uniprot_acc</a>");#Ajout de la variable (url) dans le tableau (Iteration)		
		}			
	@tab_string_url=uniq(@tab_string_url);
	$string_url = join(' ', @tab_string_url);
	print "Reconstruction terminee \n";

########################################################################
#                               PDB                                    #
########################################################################
	print "Extraction d'IDentifiants PDB....................";
	my @list_pdbID=();
		foreach my $uniprot_id (@uniprot_ID){#Parcours du tableau
			my $XML_query = qq(<?xml version="1.0" encoding="UTF-8"?><orgPdbQuery><queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
			<description>Simple query for a list of Uniprot Accession IDs:'$uniprot_id'</description><accessionIdList>'$uniprot_id'</accessionIdList></orgPdbQuery>);
		
			#----------Creation de la requête----------#
																						 
			my $request = HTTP::Request->new( POST => 'http://www.rcsb.org/pdb/rest/search/');
			$request->content_type( 'application/x-www-form-urlencoded' );
			$request->content($XML_query);

			$response = $ua->request( $request );# Retour d'un objet de type Response. Il contient differents champs 

			#----En cas d'echec de la requête-----#
				unless( $response->is_success ) { #is_success permet de savoir si la requête s'est bien deroulee
					print "\n an error occurred: ", $response->status_line, "\n";#Dans le cas contraire, affichage d'un message d'erreur avec status_line
				}

			#----En cas de reussite-----#
			my $rep=$response->content;chomp($rep); #content: contenant le texte de la page
			@list_pdbID = split("\n", $rep);
			@list_pdbID = map( { substr($_,0,4) } @list_pdbID );#Suppression des 2 derniers caractères de chaque element du tableau (ex: 1B22:1->1B22)
			$list_pdbID = join(' ', @list_pdbID);
			print "Extraction terminee \n";

########################################################################
#                            Quick GO                                  #
########################################################################
			print "Recuperation des IDentifiants Quick GO .............";
			#-----Initialisation de l'agent-----#
			my $ua = LWP::UserAgent->new;
			#-----Creation de la requête--------#
			#La méthode GET ajoute les données à l'URL. 
			my $req = HTTP::Request->new(GET => "http://www.ebi.ac.uk/QuickGO/GAnnotation?protein=".$uniprot_id."&format=tsv&col=goID,goName,aspect");
			
			#-----Execution de la requête et reception de  la reponse (stockage dans un fichier temporaire)--------#
			my $res = $ua->request($req, "temp");
			#-----Declaration de tableaux vides--------#	
			my @tab_quickgo_process=();
			my @tab_quickgo_function=();
			my @tab_quickgo_component=();
			my @tab_quickgo=();
			#----Ouverture et parcours du fichier temp--------#	
			#	my $head = <FILE>;
			open (FILE, 'temp');
				while (<FILE>) {
					chomp;
				#----Eclatement de la chaines du fichier et attribution dans une variable --------#	
					my ($goID, $goName, $aspect) = split(/\t/);
				#----Trie des Go IDs en fonction de l'Aspect puis attribution dans un tableau--------#	
						if($aspect eq "Process"){
							push(@tab_quickgo_process,"<a href=https://www.ebi.ac.uk/QuickGO/GSearch?q=$goID&what=GO title=$goName>$goID</a>");
						}
						elsif($aspect eq "Function"){
							push(@tab_quickgo_function,"<a href=https://www.ebi.ac.uk/QuickGO/GSearch?q=$goID&what=GO title=$goName>$goID</a>");
						}
						else {
							push(@tab_quickgo_component,"<a href=https://www.ebi.ac.uk/QuickGO/GSearch?q=$goID&what=GO title=$goName>$goID</a>");
						}
					#----Suppression de la redondance (package uniq) et trie (sort) des Go IDs dans l'ordre croissant--------#	
						@tab_quickgo_process = sort (uniq(@tab_quickgo_process));
						@tab_quickgo_function = sort(uniq(@tab_quickgo_function));
						@tab_quickgo_component = sort(uniq(@tab_quickgo_component));
				}
		#----Construction d'une chaine de caractère avec les elements du tableau--------#	
			$tab_quickgo_process = join(' ', @tab_quickgo_process);
			$tab_quickgo_function = join(' ', @tab_quickgo_function);
			$tab_quickgo_component = join(' ', @tab_quickgo_component);
		#----Fermeture du fichier temporaire--------#		
			close FILE;
		}
	print "Extraction terminee \n";

########################################################################
#                            NCBI                                      #
########################################################################

	my $factory = Bio::DB::EUtilities->new(-eutil  => 'esearch',
											-db    => 'gene',#Database:gene
											-term  => $gene_symbol.'[sym] AND '.$species.'[organism]',#Recherche realisee à l'aide du Gene symbol et de l'organism
											-email => 'gotama@hotmail.com',
											-retmax => 10,#Limite du nombre de resultat retourne
											-verbose=>"0");#Eviter l'affichage d'informations
	my (@id,@gi,@ids,@refseq,@refseq_p) ;  #Declaration des tableaux                                  
	@id = $factory->get_ids;#Recupère le gene ID dans un tableau
	$gene_id=join(' ', @id);
#----------------------------------------------------------------------#

	my $fetcher = Bio::DB::EUtilities->new(
										-eutil		=> 'esummary',
                                          -id      	=> \@id,
										  -db     	=> 'gene',										 
                                          -retmode 	=> 'ans1',#Détermine le type du fichier de sortie
											-email  => 'gotama@hotmail.com',
											-verbose=>"0" 
											);
	my @genename = $fetcher->get_Response->content;
	my $gene_name = join(' ', @genename);
	@genename = split(/\n/, $gene_name);
	print "Extraction du Official Full Name.....................";
	#-------------------Extraction d'elements avec grep------------#
	#L'official full name est encadre par 2 balises <NomenclatureName>.....</NomenclatureName>
	$off_full_name = first {/<NomenclatureName>/i} @genename;
	$off_full_name =~ s/<\/?NomenclatureName>//g;##Suprression des balises afin de ne conserver que leur contenu . ?present 0 ou 1 fois,g toutes les occurences
		if($off_full_name == ""){# Dans le cas où il n'y aurait pas de NomenclatureName
			$off_full_name = first { /<Description>/i } @genename;#Recherche du pattern jusqu'au premier match. Ici un grep aurait parcouru tout le fichier XML
			$off_full_name =~ s/<\/?Description>//g;#\ present ?0 ou 1 fois,g toutes les occurences; S pour Substituton
		}
	print "Extraction terminee \n";

		
#----------------------------------------------------------------------#
	my $handle = Bio::DB::EUtilities->new(
											-eutil  => 'elink',
											-email  => 'gotama@hotmail.com',
											-db     => 'nuccore',
											-dbfrom => 'gene',
											-id     => \@id, 
											-verbose=>"0",);
	@ids = $handle->get_ids;	
	my $ids1 = join(' ', @ids);
	@gi=split(/\t/,$ids1);
#----------------------------------------------------------------------#
	print "Extraction des IDentifiants Refseq - Transcrits ............";

	my $handle2 = Bio::DB::EUtilities->new(-eutil  	 => 'efetch',
												-db      => 'nuccore',
												-id      => \@gi,
												-email   => 'gotama@hotmail.com',
												-rettype => 'acc',
												-verbose=>"0");
 
	my @accs = split(m{\n},$handle2->get_Response->content);
	@accs = sort(uniq(@accs));
	my $regexpr=qr /^[NX]M/;
		foreach my $cle3 (@accs){ #Parcours du tableau
			if ($cle3 =~ $regexpr){ #SI la variable commence par NM ou XM
				push(@refseq,$cle3); #Insertion dans un tableau @refseq (Iteration)
			}
		}
	
	$refseq_mrna = join(' ', @refseq);
	print "Extraction terminee \n";
#----------------------------------------------------------------------#
	print "Extraction des IDentifiants Refseq - Proteine ............";
	
		foreach $id (@id){
		
			my $handle = Bio::DB::EUtilities->new(-eutil  => 'elink',
											  -email  => 'gotama@hotmail.com',
											  -db     => 'protein',
											  -dbfrom => 'gene',											  
												-id     => \@id,
											  -verbose=>"0");
	        ###Recuperation des GI numbers###                               										  
			@ids = $handle->get_ids;	
			my $ids2 = join(' ', @ids);
			@gi=split(/\t/,$ids2);
		###Recherche de la lsite de GI dans la base proteine###                               
			my $handle3 = Bio::DB::EUtilities->new(-eutil   => 'efetch',
										   -db      => 'protein',
                                           -id      => \@gi,
                                           -email   => 'gotama@hotmail.com',                                      
                                           -rettype => 'acc',
                                            -verbose=>"0");
			my @accs2 = split(m{\n},$handle3->get_Response->content);
			@accs2 = sort(uniq(@accs2));
			my $regexpp=qr /^[NX]P/;
				foreach my $cle4 (@accs2){
					if ($cle4 =~ $regexpp){
						push(@refseq_p,$cle4);
					}
				}
		}
	$refseq_prot = join(' ', @refseq_p);
		print "Extraction terminee \n\n";
	
#----------------------------------------------------------------------#

	print FILE_RESULT " 
		<tr>
			<td>$gene_symbol</td>
			<td>$species</td>			
			<td>$gene_id</td>
			<td>$off_full_name</td>
			<td>$ensg</td>
			<td>$genebrows_url</td> 
			<td>$xrefrna</td>
			<td>$xrefprot</td>
			<td>$refseq_mrna</td>
			<td>$refseq_prot</td>		
			<td>$id_unip</td>
			<td>$protein_name_unip</td>
			<td>$list_pdbID</td>
			<td>$string_url</td>
			<td>$tab_quickgo_component</td>
			<td>$tab_quickgo_function</td>
			<td>$tab_quickgo_process</td>
		</tr>";#Print de toutes les valeurs dans un tableau HTML
	
}
#---------------------Insertion des fichiers css, javascript, et appel de la DataTable-----------------------#
	print "Creation du tableau de sortie\n";

my $end=
	'<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.11/css/dataTables.jqueryui.min.css" />
	<link rel="stylesheet" href="https://jquery-ui.googlecode.com/svn/tags/1.8.23/themes/dark-hive/jquery-ui.css" />
	<link rel="stylesheet" href="https://cdn.datatables.net/colreorder/1.3.1/css/colReorder.dataTables.min.css" />
	<link rel="stylesheet" href="https://cdn.datatables.net/fixedheader/3.1.1/css/fixedHeader.jqueryui.min.css" />
	<link rel="stylesheet" href="https://cdn.datatables.net/responsive/2.0.2/css/responsive.jqueryui.min.css" />
	<link rel="stylesheet" href="https://cdn.datatables.net/buttons/1.1.2/css/buttons.jqueryui.min.css" />
	<link rel="stylesheet" href="http://www.thesignintelligence.com/dist/datatable/extensions/ColVis/css/dataTables.colvis.jqueryui.css" />

	<script src="https://code.jquery.com/jquery-1.12.0.min.js"></script> 
	<script src="https://cdn.datatables.net/1.10.11/js/jquery.dataTables.min.js"></script> 
	<script src="https://cdn.datatables.net/1.10.11/js/dataTables.jqueryui.min.js"></script>
	<script src="https://cdn.datatables.net/colreorder/1.3.1/js/dataTables.colReorder.min.js"></script>
	<script src="https://cdn.datatables.net/buttons/1.1.2/js/dataTables.buttons.min.js"></script>
	<script src="https://cdn.datatables.net/fixedheader/3.1.1/js/dataTables.fixedHeader.min.js"></script>
	<script src="https://cdn.datatables.net/responsive/2.0.2/js/dataTables.responsive.min.js"></script>
	<script src="https://cdn.datatables.net/responsive/2.0.2/js/responsive.jqueryui.min.js"></script>
	<script src="https://cdn.datatables.net/buttons/1.1.2/js/dataTables.buttons.min.js"></script>
	<script src="https://cdn.datatables.net/buttons/1.1.2/js/buttons.jqueryui.min.js"></script>
	<script src="http://cdn.datatables.net/buttons/1.1.2/js/buttons.print.min.js"></script>
	<script src="http://www.bacubacu.com/colresizable/js/colResizable-1.5.min.js"></script>
	<script src="http://www.thesignintelligence.com/dist/datatable/extensions/ColVis/js/dataTables.colVis.js"></script>
	<script src="https://cdn.datatables.net/buttons/1.1.2/js/buttons.html5.min.js"></script>
	<script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/2.5.0/jszip.min.js"></script>
	<style>
	.ColVis {
			text-align: center;
		}


	</style>
	<script type="text/javascript" class="init">
	$(document).ready(function(){
		var table =$(".display").DataTable({
			jQueryUI: true,
			colReorder: true,
			 fixedHeader: true,
				 responsive: true,
				 paging:false,
			 
			   

		
	});	
	new $.fn.dataTable.Buttons( table, {
		
		buttons: [
			  {
				 extend: "collection",
				 text: "Export",
				 buttons: [ "pdfHtml5", "csvHtml5", "copyHtml5", "excelHtml5" ]
			  }
		   ]
			} );
	var colvis = new $.fn.dataTable.ColVis( table );
		table.buttons( 0, null ).container().prependTo(
			table.table().container()
		);
		$( colvis.button() ).insertAfter("div.dataTables_filter");
				$(".display").colResizable({liveDrag:true});
				

	} );
	</script>


';

print FILE_RESULT "</tbody></table>";#Print fin du tableau en HTML
print FILE_RESULT "$end\n";#Print des appels des differents plugin Javascript et CSS
print FILE_RESULT"</body>
</html>";#Print fin du fichier  HTML


close FILE_RESULT; #Fermeture du fichier
