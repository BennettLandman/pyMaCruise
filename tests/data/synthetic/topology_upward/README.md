Selected Java v3.2 topology-correction references

Each NPZ contains float32 input and independent Java output in x/y/z index
order. Source AlgorithmTopologyCorrection revision 1.7 from the supplied
`toads-cruise` tree uses intensity initialization, background->object, (6,18),
regularization 0.0001 and threshold stopping disabled. Java receives no paint
mask. The existing MIPAV in-memory framework shims substitute only image and
UI storage. The unchanged selected TopologyPropagation and BinaryTree run the
numerical algorithm. Details and source hash are in provenance.json.

These references validate that configured scalar operation, not complete
workflow equivalence. Patient data and compiled Java programs are excluded.
