#include <stdio.h>
#include <stdlib.h>
#include "tinyxml2.h"

using namespace tinyxml2;

int main(int argc,char** argv){
    XMLDocument doc; 
    doc.LoadFile( "systemSetting.xml" );

    XMLElement* element= doc.FirstChildElement("systemSetting")->FirstChildElement("attribute");
    for(;element!=NULL;element=element->NextSiblingElement()){
        printf("%s:",element->FirstChildElement("title")->GetText());
        printf("%s\n",element->FirstChildElement("title")->NextSiblingElement()->GetText());
    }



    // Navigate to the title, using the convenience function,
    // with a dangerous lack of error checking.
    //const char* title = doc.FirstChildElement( "PLAY" )->FirstChildElement( "TITLE" )->GetText();
    //printf( "Name of play (1): %s\n", title );
    //
    //// Text is just another Node to TinyXML-2. The more
    //// general way to get to the XMLText:
    ////XMLText* textNode = doc.FirstChildElement( "PLAY" )->FirstChildElement( "TITLE" )->FirstChild()->ToText();
    ////title = textNode->Value();
    ////printf( "Name of play (2): %s\n", title );
    //XMLElement* element=doc.FirstChildElement("PLAY");
    //element=element->NextSiblingElement();
    //printf( "Name of play (2): %s\n", 
    //        element->FirstChildElement("TITLE")->GetText());
    return 0;
}
