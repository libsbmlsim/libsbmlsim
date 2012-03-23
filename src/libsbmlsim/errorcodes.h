#ifndef LibSBMLSim_ErrorCodes_h
#define LibSBMLSim_ErrorCodes_h

typedef enum {
  NoError,
  Unknown,
  FileNotFound,
  InvalidSBML,
  SBMLOperationFailed,
  InternalParserError,
  OutOfMemory,
  SimulationFailed
} LibsbmlsimErrorCode;

#endif  /* LibSBMLSim_ErrorCodes_h */
