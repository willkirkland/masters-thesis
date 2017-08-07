/** @file xmlValidate.c Validate XML file against schema */
/** @author Will Kirkland (UTK) */
/*  Called by xml_data_mod.f90 */

#define LIBXML_SCHEMAS_ENABLED
#include <libxml/xmlschemastypes.h>
#include <libxml/xmlschemas.h>
#include <string.h>
#include "x2nSchema.h" /* generated in Makefile from x2nSchema.xsd */ 
#include "x2nSchema_new.h" /* generated in Makefile from x2nSchema_new.xsd */ 

void xmlvalidate_(const char* xmlFile, int* length, int* returnVal)
{
	char filename[*length+1]; /*!< name of xml file          */
    xmlSchemaParserCtxtPtr contextPtr;  /*!< used by libxml2 */
    xmlSchemaPtr schemaPtr;             /*!< used by libxml2 */
    xmlSchemaValidCtxtPtr validCtxtPtr; /*!< used by libxml2 */
    int i; /*!< index variable */

    /* Convert xml filename to C string */
    for(i=0; i< *length; i++)
    {
        filename[i] = xmlFile[i];
    }
    filename[*length] = '\0';

    /* Convert xml schema filename to valid context pointer and try 
	 * validating against first schema */
    contextPtr = xmlSchemaNewMemParserCtxt(x2nSchema_xsd, x2nSchema_xsd_len);
    schemaPtr = xmlSchemaParse(contextPtr);
    validCtxtPtr = xmlSchemaNewValidCtxt(schemaPtr);
	*returnVal = xmlSchemaValidateFile(validCtxtPtr, filename,0);

	/* Return 0 if first schema validates */
	/* If first schema failes to validate, try again for second schema */
    if (*returnVal){
      contextPtr = xmlSchemaNewMemParserCtxt(x2nSchema_new_xsd, x2nSchema_new_xsd_len);
      schemaPtr = xmlSchemaParse(contextPtr);
      validCtxtPtr = xmlSchemaNewValidCtxt(schemaPtr);

	  *returnVal = xmlSchemaValidateFile(validCtxtPtr, filename, 0);
	} else return;
    
    /* Free memory used in schema validation */
	xmlSchemaFreeParserCtxt(contextPtr);
	xmlSchemaFree(schemaPtr);
	xmlSchemaFreeValidCtxt(validCtxtPtr);

	/* Return 1 if second schema validates */
	if (*returnVal){
	  *returnVal = -1;
	  return;
	} else{
	  *returnVal = 1;
	  return;
	}
}
