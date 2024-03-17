#set page(
  numbering: "1"
)

// #show ref: set text(fill: rgb(0, 0, 255))

#set page(margin: 1.75in)
#set par(leading: 0.55em, justify: true)
#set text(font: "New Computer Modern")
#show raw: set text(font: "New Computer Modern Mono")
#show heading: set block(above: 1.4em, below: 1em)

= Structural biology validation landscape

Structural biology, as a branch of molecular biology, is a life science that has
seen much impact in the scientific world. This is proven by the fact that as of
now (March 2024) there were 81 Nobel prizes awarded in relation to it
@biology-nobel. It has contributed with methods for examining and
acquiring biomacromolecular scientific data. These data constitute structures
(i.e. coordinate models
@coord-models of polymers that are biological in origin (i.e. nucleic
acids, polysacharrides, and proteins).

These assembled models in computer-readable form are frequently stored in
databases which aggregate results obtained by various experimental methods. The
largest and most notable curated database is the Protein Data Bank (PDB) @pdb.  Currently, this freely-accessible archive stores over 217,000 entries
@pdb-entry-stats.  Most data in the PDB is created using methods such as Nuclear
Magnetic Resonance (NMR) spectroscopy, X-ray crystallography, or Electron
Microscopy (EM) @pdb-methods.

Because of the experimental nature of these methods, this data can contain
errors caused by the imperfections of measurements or simply, human error
@structure-errors[p. 1]. That is why significant effort has been made to develop
a consensus on how to validate this data. The PDB consortium has organised
so-called Validation Task Forces (VTFs), i.e. groups of experts on each
experimental method used in developing structure models @pdb-validation[p. 1917]. Based on the results of
the task forces, the PDB started to offer validation reports for most of the
structures in the database @pdb-validation[p. 1917]. These reports are generated by the Validation
Pipeline. The pipeline modules run validations on the structure using existing
bioinformatical software @pdb-validation[p. 1921]. After running the validations, the results are
aggregated and the reports are composed into a human-readable form (PDF) and a
computer-readable form (XML) @pdb-validation[p. 1921]. When a new structure is deposited into the PDB
database, the pipeline is run and the reports are stored alongside the deposited
structure @pdb-validation[p. 1916].

To allow users to check their experimental data without beginning the deposition
process, the PDB consortium created the wwPDB standalone validation web server @pdb-validation[p. 1917].
It is a scaled-down version of the validation pipeline that is used during
deposition proper @pdb-standalone-server.

// TODO: Mention the 600 weekly uses of the API in 2017
Public access to the standalone server is provided using a web server, a CLI
application, and a Python library @pdb-validation[p. 1920] @pdb-standalone-server-details. However, the service is constituted of only 2
compute servers (2 x 16 CPUs, 2 x 128 GB RAM) for the entire world. The low
throughput make it impractical to be used in research projects that need to
validate a lot of structures. Even though the software components of the
validation pipeline are public @pdb-validation[p. 1922], the implementation of the pipeline itself is
not. For this reason, it is not possible to deploy it using on-premise
infrastructure.

To allow for simple deployment and high-throughput validations, a new solution
is needed. Unfortunately, not all components of the pipeline are publicly
available. In this work, we reimplement the publicly available parts of the
validation pipeline using software that allows for easy scalability and
deployment.

Since the PDB validation pipeline uses multiple software components with various
dependencies, it is necessary to create an environment in which all of the
components can run. To achieve this goal, we use containers, which guarantee
that the software suite can run in any cloud or non-cloud environment
@containers[p. 2-3]. For packaging the components and their dependencies, we use
_Docker_, a set of tools that allows easy creation of containers @containers[p. 3].

To allow for scalability, _Kubernetes_ is utilized. An open-source system
originally developed at _Google_ allows for automatic deployment and scaling of
containerized applications @kubernetes.

To simplify deployment, we use _Ansible_, an open-source automation tool
@ansible which allows for deployment in any _Kubernetes_ cluster.

#bibliography("bibliography.yaml", style: "ieee")