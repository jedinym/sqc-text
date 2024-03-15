#set page(
  numbering: "1"
)

= Structural biology validation landscape

Structural biology, as a branch of molecular biology, is a life science that has
seen much impact in the scientific world. This is proven by the fact that as of
now (March 2024) there were 81 Nobel prizes awarded in relation to it
#cite(<biology-nobel>). It has contributed with methods for examining and
acquiring biomacromolecular scientific data. These data constitute structures
(i.e. coordinate models
#cite(<coord-models>)) of polymers that are biological in origin (i.e. nucleic
acids, polysacharrides, and proteins).

These assembled models in computer-readable form are frequently stored in
databases which aggregate results obtained by various experimental methods. The
largest and most notable curated database is the Protein Data Bank (PDB).
Currently, this freely-accessible archive is storing over 217,000 entries.

Most data in the PDB is created using methods
such as Nuclear Magnetic Resonance (NMR) spectroscopy, X-ray crystallography, or
Electron Microscopy (EM). Because of the experimental nature of these methods,
this data can contain errors caused by the imperfections of measurements or
simply, human error. That is why significant effort has been made to develop a
consensus on how to validate this data. The PDB consortium has organised
so-called Validation Task Forces (VTFs), i.e. groups of experts on each
experimental method used in developing structure models. Based on the results of
the task forces, the PDB started to offer validation reports for most of the
structures in the database.

To allow users to check their experimental data without beginning the deposition
process, the PDB consortium created the wwPDB validation service. It is
a scaled-down version of the validation pipeline that is used during deposition
proper.

Access to the service is provided using a web server, that is intended for human
interaction. Furthermore,the service is constituted of only 2 compute servers (2
x 16 CPUs, 2 x 128 GB RAM) for the entire world. The lack of simple programmatic
access and low throughput make it impractical to be used in research projects
that need to validate a lot of structures.

#bibliography("bibliography.yaml")