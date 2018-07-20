# How to build a database
You need three kind of files:

## fasta sequences (with the extension .fas) 
They can easily be found if you get the right uniprot ID with wget http://www.uniprot.org/uniprot/[UNIPROT_ID].fasta
You can also create them reading the .pdb structure file of the protein.
Put these files in the folder [project_folder]/dataset/fas/

## Mutiple Sequence Alignments (with the extension .aln)
To create them: use hh-suite (documentation online)
hhblits -i [your_fasta_file] -opsi [your_psi_file] -cpu 4 -n 1 -d [june 2013 uniprot database] -d [june 2016 uniprot database]
A lot of them should be available. Their computation take ~3mn per sequence on graham.

Then reformat them in the .sto format using the reformat.pl script from hhsuite as well:

for file in `ls *.psi`; do
reformat.pl $file `echo $file | awk -F . '{print $1}'`.sto
done

This step can be long. You can delete the .psi files after or move them elsewhere.
The format then is not yet really practical for usage. The little routine reformat_sto.sh that can be found in the folder [project_folder]/dataset/routines/ solves the issue. Put all you sto files in the folder [project_folder]/dataset/aln/ and then type:
./reformat_sto.sh [project_folder]/dataset/
make sure to type the whole path (not . ) and to put a / after dataset. It's very important. This step can be pretty long.
Now the files should be ready for the final step and in the folder [project_folder]/dataset/aln/

## structure files (most likely compressed with the extension pdb.xz)
If you need ModBase files, you have to send a mail to modbase@salilab.org with the list of models you need.
To extract the contact maps from these structure files, you should put all your structures in the folder [project_folder]/dataset/pdbxz, go into the folder [project_folder]/scripts, open the file setup.py and set the variable folder to [project_folder] then launch the script create_cmap.py

Multiple errors can occur due to the lack of a unity between the structure of PDB files in themselves. For instance if the structure possess chains you might want to replace line 19 and 20 by:
                if line.split()[0] == 'ATOM' and line.split()[2]=='CA' and line.split()[5]==[letter of the chain you want]:
                    _x, _y, _z = line.split()[6:9]
Some files are broken in some way. You will always end up with files with mistakes. They'll usually wont cause any trouble to the algorithm and you can find their list in the file [project_folder]/database/log/badfiles.log

## Finally, launch [project_folder]/scripts/serializer.py to finish the database building
