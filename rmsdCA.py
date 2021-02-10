# load library
from pymol import cmd

# main function
@cmd.extend
def rmsdCA(refMol, tgtMol, puttyview="Y", bfactorca="CA"):
    """
    Written up by Yufeng Tong 2020-01-29
    Inspired by rmsdByRes Zhenting Gao on 7/28/2016
    and by loadBfacts from PyMOLWiki

    USAGE

      rmsdCA refMol, tgtMol, puttyview="Y", bfactorca="CA"

      Calculate the RMSD of Calpha and backbone for each residue pairs from
        two chains of proteins with same residue numbers and close
        conformations for the two structures. The two proteins can be 
        mutants. Only residue numbers but not identity are matched.
      It will generate a putty representation of the target molecule by 
        default. You can specify whether to use the RMSD of Calpha or Backbone
        for the coloring of the putty display.
      I think it's not neccessary to calculate RMSD of the whole residue.
      A CSV file is saved of the calculated RMSD of Calpha and backbone atoms.
      A PDB file is saved with BFactor replaced by the RMSD.

    Workflow
      Read reference and target pdb files
      1. Select two objects for comparison
      2. Align two structures
      3. rmsdByRes refMol, tgtMol
    
    TODO:
      1. Save output files to the folder pdb files are loaded.
      2. Handle missing residues.
    """


# Create temporary objects, exclude alternative conformation B
    cmd.delete("tgt_gzt rmsdPuttyScale")
    cmd.create("ref_gzt", refMol + " and polymer and not alt B")
    cmd.alter("ref_gzt", "chain='A'")
    cmd.alter("ref_gzt", "segi=''")
    cmd.create("tgt_gzt", tgtMol + " and polymer and not alt B")
    cmd.alter("tgt_gzt", "chain='A'")
    cmd.alter("tgt_gzt", "segi=''")
    cmd.alter("tgt_gzt", "b=-1.0")

# parameters
    outputText = ""
    csvHeadline = ("refMol,refRes,refResID,target,residueName,residueId,"
                   "rmsdResCa,rmsdResBackbone\n")
    outputFile = 'rmsd%s_%s.csv' % (bfactorca, tgtMol)
    outputpdb = 'rmsdBFactor%s_%s.pdb' % (bfactorca, tgtMol)
    bfactors = []

# select alpha carbon of selected residues in reference structure
    calpha_ref = cmd.get_model("ref_gzt and name CA")
    calpha_tgt = cmd.get_model("tgt_gzt and name CA")

# build an index dict of {resi : index}
    ref_dict = {int(v.resi): k for k, v in enumerate(calpha_ref.atom)}
    tgt_dict = {int(v.resi): k for k, v in enumerate(calpha_tgt.atom)}

# loop through residues in the target molecule and calculate per residue RMSDs.
    for k in sorted(ref_dict):
        if k in tgt_dict:
            ref_atom = calpha_ref.atom[ref_dict[k]]
            tgt_atom = calpha_tgt.atom[tgt_dict[k]]
            rmsdResCA = cmd.rms_cur(
                "ref_gzt and n. ca and i. " + ref_atom.resi,
                "tgt_gzt and n. ca and i. " + tgt_atom.resi, matchmaker=-1)
            rmsdResBb = cmd.rms_cur(
                "ref_gzt and n. ca+n+c+o and i. " +
                ref_atom.resi,
                "tgt_gzt and n. ca+n+c+o and i. " +
                tgt_atom.resi,
                matchmaker=-
                1)
            if bfactorca=="bb":
                bfactor = rmsdResBb
            else:
                bfactor = rmsdResCA
            cmd.alter("%s and i. %s" % ("tgt_gzt", tgt_atom.resi),
                      "b=%s" % bfactor)
            bfactors.append(bfactor)

            outputText += "%s,%s,%s,%s,%s,%s,%.3f,%.3f\n" % (
                refMol, ref_atom.resn, ref_atom.resi, tgtMol, tgt_atom.resn,
                tgt_atom.resi, rmsdResCA, rmsdResBb)

    print(outputText)
    max_b, min_b = max(bfactors), min(bfactors)
    print("RMSD: max:%.3f; min:%.3f" % (max_b, min_b))

    if puttyview == "Y":
        cmd.show_as("cartoon", "tgt_gzt")
        cmd.cartoon("putty", "tgt_gzt")
        cmd.set("cartoon_putty_scale_min", min_b,  "tgt_gzt")
        cmd.set("cartoon_putty_scale_max", max_b, "tgt_gzt")
        cmd.set("cartoon_putty_transform", 0, "tgt_gzt")
        cmd.spectrum("b", "rainbow", "tgt_gzt and n. CA")
        cmd.ramp_new(
            "rmsdPuttyScale", "tgt_gzt", [min_b, max_b], "rainbow")
        cmd.recolor()

# Destroy temporary objects
    cmd.delete("ref_gzt")

# Save data into csv and pdb
    f = open(outputFile, 'w+')
    f.write(csvHeadline)
    f.write(outputText)
    f.close()
    cmd.save(outputpdb, "tgt_gzt")

    print("Results are saved in %s and %s " % (outputFile, outputpdb))
