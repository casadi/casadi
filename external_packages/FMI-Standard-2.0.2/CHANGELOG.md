# Changelog

## v2.0.2 (2020-12-15)

*	`p. 22` Fixed: Discourage use of memory management functions (#931)
*	`p. 27` Fixed formula for directional derivatives (#668)
*	`p. 29` Fixed: Correct code example describing Jacobian construction using fmi2GetDirectionalDerivative (#867)
*	`p. 29` Fixed: Formatting of hints on graph coloring algorithm  (#703)
*	`p. 37ff` Fixed: Unit speficiation error relative quantity usage (#1267)
*	`p. 67` Fixed: "Apostrophe" instead of "Right Single Quotation Mark" in EBNF (#692)
*	`p. 68` Clarified: allowing 0-based arrays (#1197)
*	`p. 70` Added: /extra directory as clarification of storage location for additional files, non-functional change (#1207)
*	`p. 91` Fixed: typos in pseudocode (#1269)

## v2.0.1 (2019-10-31)

Specification document:

* improved wording in the whole document
* typographical errors fixed in the whole document (#444, #623)
* `p. 3` Clarified BSD version in standard documents (#347, #494)
* `p. 8` Remove references to outdated ImplementationHint document (#352)
* `p. 14` clarified:  Who shall define the prefix for source code FMUs? (#420)
* `p. 27 ff., p. 88 ff.` fixed pseudo code examples:
  * fmi2GetDirectionalDerivative: Inconsistency in specification and pseudocode example (#311)
  * detect a state change during event iteration: wrong in pseudo code  (#333)
  *  Wrong call sequence for fmi2SetTime in ME pseudo code example (#296)
  *  Fix Pseudocode Example (#308)
* `p. 31` Improved ModelStructure description in FMI 2.0 specification (#374)
* `p. 22` fixed: inconsistency between description and state machine (#338)
* `p. 33, p.68` provide license information (#417)
* `p. 33` Clarifying FMI version string: forward/backwards compatibility; version string (#388, #498)
* `p. 41` Clarified: do not use empty unit definitions (#337)
* `p. 41` Clarified that all units must be defined in UnitDefinitions (#337, #500)
* `p. 49` Remove reference to fmi2ExitEventMode  (#432, #450)
* `p. 52 ff.` Clarification of alias variables (and variability) (#436)
* `p. 55` clarified: Changing parameter values in the modeldescription.xml file discussion specification (#387)
* `p. 60`  Clarified the definition of the number of continuous-time state variables model-exchange (#396)
* `p. 65` added example for Dependency information for FMUs for ModelExchange and CoSimulation (#380)
* `p. 66` clarification:  Setting input start values before fmi2EnterInitializationMode (clarification in 2.0.1) (#440)
* `p. 66` Clarification that variable names must be unique  (#415)
* `p. 67` Clarified:  Are all of the API functions mandatory? (#411)
* `p. 68` Update text under state chart in FMI 2 (#431)
* `p. 68` clarified: Backslash as path separator not be allowed in FMUs (#547)
* `p. 68` clarified:  Handling international encodings (#95)
* `p. 68` clarified: Location of additional DLLs (#401)
* `p. 69` clarified: Order of scalar variables with respect to structured naming (#413)
* `p. 73` Improve sentence using nextMode (#393, #456)
* `p. 77` fixed: Inconsistent definition of the calling conditions for fmi2SetContinuousStates (#407)
* `p. 80` Clarification for pure discrete FMUs
  * Pure discrete-time FMU with state events (#405)
  * Non-determinism regarding the pure discrete-time FMU variant (#409)
* `p. 82` Fixed:  Description of fmi2GetContinuousStates is out of date (#334, #501)
* `p. 85` inconsistent definitions of when fmi2GetDerivatives and fmi2GetNominalsOfContinuousStates may be called (#301)
* `p. 87` fixed: fmi2SetTime in Event Mode of Model Exchange (#394)
* `p. 87` fixed: Can fmi2SetReal be used for changing the value of continuous-time state variables? (#391, #392, #402)
* `p. 99` Clarified mathematical model for co-simulation (#460)
* `p. 103` fmi2CancelStep clarification (#426, #460)
* `p. 115` Fixed literature links (#478)

C headers (no functional changes):

* Add parameter names to function prototypes (#2)
* License changed to 2-clause BSD (w/o extensions) (#347, #494)
* fixed: fmiCallbackFunction structure declare const function pointers (#216)

XSD schema (no functional changes)

* Clarifying FMI version string in xsd documentation string (#388 , #498)

Distribution archive:

* file names and structure changed
