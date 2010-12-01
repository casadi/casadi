#ifndef DUMP_HPP
#define DUMP_HPP

#include "tinyxml.h"

const char * getIndent( unsigned int numIndents );

/** \brief  same as getIndent but no "+" at the end */
const char * getIndentAlt( unsigned int numIndents );

int dump_attribs_to_stdout(TiXmlElement* pElement, unsigned int indent);

void dump_to_stdout( TiXmlNode* pParent, unsigned int indent = 0 );

/** \brief  load the named file and dump its structure to STDOUT */
void dump_to_stdout(const char* pFilename);


#endif