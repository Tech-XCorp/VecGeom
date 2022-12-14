//!    \file Helper.cpp
//!    \brief Defines helper functions to manipulate data from XML
//!
//!    \authors Author:  Dmitry Savin <sd57@protonmail.ch>
//!
//!    \license Distributed under the Apache license 2.0

#include "Helper.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <array>
#include <algorithm>

#include "xercesc/dom/DOMNamedNodeMap.hpp"
#include "xercesc/dom/DOMNode.hpp"
#include "xercesc/dom/DOMAttr.hpp"
#include "xercesc/dom/DOMTypeInfo.hpp"
#include "xercesc/dom/DOMElement.hpp"

namespace {
static const std::array<std::string, 13> nodeTypeNames{
    {"invalid", // the numbering starts from 1
     "ELEMENT_NODE", "ATTRIBUTE_NODE", "TEXT_NODE", "CDATA_SECTION_NODE", "ENTITY_REFERENCE_NODE", "ENTITY_NODE",
     "PROCESSING_INSTRUCTION_NODE", "COMMENT_NODE", "DOCUMENT_NODE", "DOCUMENT_TYPE_NODE", "DOCUMENT_FRAGMENT_NODE",
     "NOTATION_NODE"}};
}

namespace vgdml {
namespace Helper {
template <>
std::string Transcode(const XMLCh *const anXMLstring);
template <>
double Transcode(const XMLCh *const anXMLstring);
template <>
int Transcode(const XMLCh *const anXMLstring);

template <>
std::string GetAttribute(std::string const &attrName,
                         XERCES_CPP_NAMESPACE_QUALIFIER DOMNamedNodeMap const *theAttributes);
template <>
double GetAttribute(std::string const &attrName, XERCES_CPP_NAMESPACE_QUALIFIER DOMNamedNodeMap const *theAttributes);
template <>
int GetAttribute(std::string const &attrName, XERCES_CPP_NAMESPACE_QUALIFIER DOMNamedNodeMap const *theAttributes);

template <>
std::string GetAttribute(std::string const &attrName,
                         XERCES_CPP_NAMESPACE_QUALIFIER DOMNamedNodeMap const *theAttributes)
{
  auto *attrXMLName = XERCES_CPP_NAMESPACE_QUALIFIER XMLString::transcode(attrName.c_str());
  auto *attribute   = theAttributes->getNamedItem(attrXMLName);
  XERCES_CPP_NAMESPACE_QUALIFIER XMLString::release(&attrXMLName);
  return attribute ? vgdml::Helper::Transcode(attribute->getNodeValue()) : "";
};

template <>
double GetAttribute(std::string const &attrName, XERCES_CPP_NAMESPACE_QUALIFIER DOMNamedNodeMap const *theAttributes)
{
  auto const strAttribure = GetAttribute(attrName, theAttributes);
  try {
    return std::stod(strAttribure);

  } catch (std::invalid_argument &) {
    return std::numeric_limits<double>::quiet_NaN();
  }
};

template <>
int GetAttribute(std::string const &attrName, XERCES_CPP_NAMESPACE_QUALIFIER DOMNamedNodeMap const *theAttributes)
{
  auto const strAttribure = GetAttribute(attrName, theAttributes);
  if (strAttribure.empty()) return 0;
  return std::stoi(strAttribure);
};

std::map<std::string const, std::string const> GetAttributes(XERCES_CPP_NAMESPACE_QUALIFIER DOMNode const *aDOMNode)
{
  auto const *const theAttributes = aDOMNode->getAttributes(); // TODO check nullptr
  return GetAttributes(theAttributes);
}

std::map<std::string const, std::string const> GetAttributes(
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNamedNodeMap const *theAttributes)
{
  std::map<std::string const, std::string const> result;
  auto const nAttributes = theAttributes->getLength();
  for (auto ind = 0u; ind < nAttributes; ++ind) {
    auto const *const anAttribute = theAttributes->item(ind);
    auto const attributeName      = Helper::Transcode(anAttribute->getNodeName());
    auto const attributeValue     = GetAttribute(attributeName, theAttributes);
    result.insert(std::make_pair(attributeName, attributeValue));
  }
  return result;
}

template <>
std::string Transcode(const XMLCh *const anXMLstring)
{
  // TODO use u16string and then c++ standard codecvt
  auto *aCstring     = XERCES_CPP_NAMESPACE_QUALIFIER XMLString::transcode(anXMLstring);
  auto const aString = aCstring ? std::string(aCstring) : "";
  XERCES_CPP_NAMESPACE_QUALIFIER XMLString::release(&aCstring);
  return aString;
}

// TODO other types, use constexpr_if or SFINAE
template <>
double Transcode(const XMLCh *const anXMLstring)
{
  auto const aString = Transcode(anXMLstring);
  auto const aDouble = std::stod(aString);
  return aDouble;
}

template <>
int Transcode(const XMLCh *const anXMLstring)
{
  auto const aString = Transcode(anXMLstring);
  auto const anInt   = std::stoi(aString);
  return anInt;
}

std::string GetNodeInformation(XERCES_CPP_NAMESPACE_QUALIFIER DOMNode const *aDOMNode)
{
  std::stringstream aStream;
  auto const theNodeName         = Transcode(aDOMNode->getNodeName());
  auto const theNodeText         = Transcode(aDOMNode->getTextContent());
  auto const *const asDOMElement = dynamic_cast<XERCES_CPP_NAMESPACE_QUALIFIER DOMElement const *>(aDOMNode);
  auto const theNodeTypeID       = aDOMNode->getNodeType();
  aStream << "node type is \"" << nodeTypeNames.at(theNodeTypeID) << "\"";
  aStream << ", node name is \"" << theNodeName << "\"";
  aStream << ", node text content is ";
  if (IsWhitespace(theNodeText))
    aStream << "whitespace";
  else
    aStream << "\"" << theNodeText << "\"";
  if (asDOMElement) {
    auto const theNodeLocalName        = Transcode(aDOMNode->getLocalName());
    auto const theChildTagName         = Transcode(asDOMElement->getTagName());
    auto const *const attributes       = aDOMNode->getAttributes();
    auto const *const typeInfo         = asDOMElement->getSchemaTypeInfo();
    auto const nodeSchemaTypeName      = Transcode(typeInfo->getTypeName());
    auto const nodeSchemaTypeNamespace = Transcode(typeInfo->getTypeNamespace());
    aStream << ", local node name is \"" << theNodeLocalName << "\"";
    aStream << ", node tag name is \"" << theChildTagName << "\"";
    aStream << ", node type name is \"" << nodeSchemaTypeName << "\"";
    aStream << ", node type namespace is \"" << nodeSchemaTypeNamespace << "\"";
    auto const nAttributes = attributes->getLength();
    aStream << ", it has " << nAttributes << " attributes ( ";
    auto const attributeMap = GetAttributes(asDOMElement);
    for (auto const &it : attributeMap) {
      auto const attributeName  = it.first;
      auto const attributeValue = it.second;
      aStream << "\"" << attributeName << "\":\"" << attributeValue << "\" ";
    }
    aStream << ")";
  }
  return aStream.str();
}

bool IsWhitespace(std::string const &aString)
{
  return std::find_if(aString.begin(), aString.end(), [](char c) { return std::isalnum(c); }) == aString.end();
}

} // namespace Helper
} // namespace vgdml
