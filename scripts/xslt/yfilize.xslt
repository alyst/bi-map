<?xml version="1.0"?>
<xsl:stylesheet id="yfilize" version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:g="http://graphml.graphdrawing.org/xmlns"
    xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:y="http://www.yworks.com/xml/graphml"
    xmlns:math="http://exslt.org/math" xmlns:str="http://exslt.org/strings"
    extension-element-prefixes="str math" exclude-result-prefixes="g">
    <!--<xsl:namespace-alias stylesheet-prefix="g" result-prefix=""/> -->
    <xsl:output method="xml" indent="no"/>
    <xsl:param name="NodeAnnotationDisplay">no</xsl:param>
    <xsl:param name="NodeStabilityDisplay">no</xsl:param>
    <xsl:attribute-set name="id">
        <xsl:attribute name="id">
                <xsl:value-of select='@id'/>
            </xsl:attribute>
    </xsl:attribute-set>
    <xsl:attribute-set name="edge_nodes">
        <xsl:attribute name="source">
                <xsl:value-of select='@source'/>
            </xsl:attribute>
        <xsl:attribute name="target">
                <xsl:value-of select='@target'/>
            </xsl:attribute>
    </xsl:attribute-set>
    <xsl:variable name='NodeLabelKey'
        select="/g:graphml/g:key[@for='node' and @attr.name='short_id']/@id"/>
    <xsl:variable name='NodeExpDescrKey'
        select="/g:graphml/g:key[@for='node' and @attr.name='experimental_description']/@id"/>
    <xsl:variable name='NodeTypeKey'
        select="/g:graphml/g:key[@for='node' and @attr.name='node_type']/@id"/>
    <xsl:variable name='NodeAvgFreqKey'
        select="/g:graphml/g:key[@for='node' and @attr.name='avg_pairs_freq']/@id"/>
    <xsl:variable name='NodeAnnotationKey'
        select="/g:graphml/g:key[@for='node' and @attr.name='annotation_type']/@id"/>
    <xsl:variable name='EdgeTypeKey'
        select="/g:graphml/g:key[@for='edge' and @attr.name='type']/@id"/>
    <xsl:variable name='EdgeHexIntensityKey'
        select="/g:graphml/g:key[@for='edge' and @attr.name='hex.intensity']/@id"/>
    <xsl:template match="/g:graphml">
        <xsl:copy>
            <xsl:value-of select='$NodeLabelKey'/>
            <xsl:copy-of select='g:key'/>
            <xsl:comment>Attributes specific for yFiles</xsl:comment>
            <key for="graphml" id="y:res" yfiles.type="resources"/>
            <key for="port" id="y:pg" yfiles.type="portgraphics"/>
            <key for="port" id="y:pg" yfiles.type="portgeometry"/>
            <key for="port" id="y:pud" yfiles.type="portuserdata"/>
            <key for="node" id="y:ng" yfiles.type="nodegraphics"/>
            <key for="edge" id="y:eg" yfiles.type="edgegraphics"/>
            <xsl:apply-templates/>
        </xsl:copy>
    </xsl:template>
    <!-- graph template -->
    <xsl:template match="g:graph">
        <xsl:copy use-attribute-sets="id">
            <xsl:apply-templates/>
        </xsl:copy>
    </xsl:template>
    <!-- node template -->
    <xsl:template match="g:node">
        <xsl:variable name='NodeLabel' select="g:data[@key=$NodeLabelKey]/text()"/>
        <xsl:variable name='NodeType' select="g:data[@key=$NodeTypeKey]/text()"/>
        <xsl:variable name='NodeFreq' select="number(g:data[@key=$NodeAvgFreqKey]/text())"/>
        <xsl:variable name='NodeAnnotation' select="g:data[@key=$NodeAnnotationKey]/text()"/>
        <xsl:variable name='NodeColor'>
            <xsl:choose>
                <xsl:when test="$NodeType='bait'">#FF2222</xsl:when>
                <xsl:when test="$NodeType='baits_cluster'">#FFAAAA</xsl:when>
                <xsl:when test="$NodeType='prey'">#55CCFF</xsl:when>
                <xsl:when test="$NodeType='proteins_cluster'">#BBCCFF</xsl:when>
                <xsl:when test="$NodeType='samples_cluster'">#FF8855</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='NodeShape'>
            <xsl:choose>
                <xsl:when test="$NodeType='bait'">ellipse</xsl:when>
                <xsl:when test="$NodeType='baits_cluster'">rectangle</xsl:when>
                <xsl:when test="$NodeType='prey'">diamond</xsl:when>
                <xsl:when test="$NodeType='proteins_cluster'">rectangle</xsl:when>
                <xsl:when test="$NodeType='samples_cluster'">rectangle</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='NodeSize'>
            <xsl:choose>
                <xsl:when test="$NodeType='bait'">50</xsl:when>
                <xsl:when test="$NodeType='baits_cluster'">60</xsl:when>
                <xsl:when test="$NodeType='prey'">30</xsl:when>
                <xsl:when test="$NodeType='proteins_cluster'">60</xsl:when>
                <xsl:when test="$NodeType='samples_cluster'">40</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='NodeFontSize'>
            <xsl:choose>
                <xsl:when test="$NodeType='bait'">15</xsl:when>
                <xsl:when test="$NodeType='baits_cluster'">18</xsl:when>
                <xsl:when test="$NodeType='prey'">9</xsl:when>
                <xsl:when test="$NodeType='proteins_cluster'">18</xsl:when>
                <xsl:when test="$NodeType='samples_cluster'">15</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='NodeLabelBackgroundColor'>
            <xsl:choose>
                <xsl:when test="$NodeType='baits_cluster'">#FFCC99</xsl:when>
                <xsl:when test="$NodeType='proteins_cluster'">#99CCFF</xsl:when>
                <xsl:when test="$NodeType='samples_cluster'">#FFCC55</xsl:when>
                <xsl:otherwise>#000000</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='NodeBorderStyle'>
           <xsl:choose>
               <xsl:when test="$NodeStabilityDisplay = 'yes'">
                <xsl:choose>
                    <xsl:when test="$NodeFreq &lt;= 0.4">dotted</xsl:when>
                    <xsl:when test="$NodeFreq &gt; 0.4 and $NodeFreq &lt;= 0.6">dashed</xsl:when>
                    <xsl:when test="$NodeFreq &gt; 0.6">line</xsl:when>
                </xsl:choose>
               </xsl:when>
               <xsl:otherwise>line</xsl:otherwise>
           </xsl:choose>
        </xsl:variable>
        <xsl:variable name='NodeBorderWidth'>
            <xsl:choose>
                <xsl:when test="$NodeStabilityDisplay = 'yes'">
                    <xsl:choose>
                        <xsl:when test="$NodeFreq &gt; 0.9">3.0</xsl:when>
                        <xsl:when test="$NodeFreq &gt; 0.7 and $NodeFreq &lt;= 0.8">2.0</xsl:when>
                        <xsl:otherwise>1.0</xsl:otherwise>
                    </xsl:choose>
                </xsl:when>
                <xsl:when test="$NodeAnnotationDisplay = 'yes' and $NodeAnnotation != ''">3.0</xsl:when>
                <xsl:otherwise>1.0</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='NodeBorderColor'>
            <xsl:choose>
                <xsl:when test="$NodeAnnotationDisplay = 'yes'">
                    <xsl:choose>
                        <xsl:when test="$NodeAnnotation = 'known'">#66FF66</xsl:when>
                        <xsl:when test="$NodeAnnotation = 'putative'">#FFFF33</xsl:when>
                        <xsl:when test="$NodeAnnotation = 'isozymes'">#FF9933</xsl:when>
                        <xsl:when test="$NodeAnnotation = 'incorrect'">#FF3333</xsl:when>
                        <xsl:otherwise>#333333</xsl:otherwise>
                    </xsl:choose>
                </xsl:when>
                <xsl:otherwise>#666699</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>
        <xsl:element name="node" use-attribute-sets="id">
            <!-- simple node style -->
            <xsl:copy-of select='g:data'/>
            <xsl:if test="not(g:graph)">
                <data key="y:ng">
                    <y:ShapeNode>
                        <y:Geometry height="{$NodeSize}" width="{$NodeSize}"/>
                        <y:Fill color="{$NodeColor}" transparent="false"/>
                        <y:BorderStyle color="#000000" type="line"
                            width="1.0"/>
                        <y:NodeLabel alignment="center"
                            autoSizePolicy="content" fontFamily="Dialog"
                            fontSize="{$NodeFontSize}" fontStyle="plain"
                            hasBackgroundColor="false" hasLineColor="false"
                            hasText="{boolean($NodeLabel)}" height="4.0"
                            modelName="internal" modelPosition="c" textColor="#000000"
                            visible="true" width="4.0" x="13.0" y="13.0">
                            <xsl:value-of select="$NodeLabel"/>
                        </y:NodeLabel>
                        <y:Shape type="{$NodeShape}"/>
                    </y:ShapeNode>
                </data>
            </xsl:if>
            <!-- group node style -->
            <xsl:if test="g:graph">
                <data key="y:ng">
                    <y:ProxyAutoBoundsNode>
                        <y:Realizers active="0">
                            <y:GroupNode>
                                <y:Fill color="{$NodeColor}" transparent="false"/>
                                <y:BorderStyle color="{$NodeBorderColor}"
                                    type="{$NodeBorderStyle}" width="{$NodeBorderWidth}"/>
                                <y:NodeLabel alignment="center"
                                    autoSizePolicy="content" borderDistance="0.0"
                                    fontFamily="Dialog" fontSize="18" fontStyle="italic"
                                    hasLineColor="false" hasBackgroundColor="false"
                                    modelName="internal" modelPosition="t"
                                    textColor="#000000" visible="true">
                                    <xsl:value-of select='$NodeLabel'/>
                                </y:NodeLabel>
                                <y:Shape type="ellipse"/>
                                <y:State closed="false" innerGraphDisplayEnabled="false"/>
                                <y:Insets bottom="0" left="0" right="0" top="0"/>
                                <y:BorderInsets bottom="0" left="0" right="0" top="0"/>
                            </y:GroupNode>
                            <y:GroupNode>
                                <y:Fill color="{$NodeColor}" transparent="false"/>
                                <y:BorderStyle color="{$NodeBorderColor}"
                                    type="{$NodeBorderStyle}" width="{$NodeBorderWidth}"/>
                                <y:NodeLabel alignment="center"
                                    autoSizePolicy="content" borderDistance="0.0"
                                    fontFamily="Dialog" fontSize="18" fontStyle="italic"
                                    hasLineColor="false" hasBackgroundColor="false"
                                    modelName="internal" modelPosition="t"
                                    textColor="#000000" visible="true">
                                    <xsl:value-of select='$NodeLabel'/>
                                </y:NodeLabel>
                                <y:Shape type="ellipse"/>
                                <y:State closed="true" innerGraphDisplayEnabled="true"/>
                                <y:Insets bottom="0" left="0" right="0" top="0"/>
                                <y:BorderInsets bottom="0" left="0" right="0" top="0"/>
                            </y:GroupNode>
                        </y:Realizers>
                    </y:ProxyAutoBoundsNode>
                </data>
                <xsl:apply-templates select="g:graph"/>
            </xsl:if>
        </xsl:element>
    </xsl:template>
    <!-- edge template -->
    <xsl:template match="g:edge">
        <xsl:variable name='EdgeType' select="g:data[@key=$EdgeTypeKey]/text()"/>
        <xsl:variable name='HexIntensity'
            select="g:data[@key=$EdgeHexIntensityKey]/text()"/>
        <xsl:variable name='LineStyle'>
            <xsl:choose>
                <xsl:when test="$EdgeType='bait_prey'">line</xsl:when>
                <xsl:when test="$EdgeType='bimap_block'">line</xsl:when>
                <xsl:when test="$EdgeType='samples_cluster_composition'">dashed</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='IsDirected'>
            <xsl:choose>
                <xsl:when test="$EdgeType='bait_prey'">true</xsl:when>
                <xsl:when test="$EdgeType='bimap_block'">true</xsl:when>
                <xsl:when test="$EdgeType='samples_cluster_composition'">false</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='Width'>
            <xsl:choose>
                <xsl:when test="$EdgeType='bait_prey'">3.0</xsl:when>
                <xsl:when test="$EdgeType='bimap_block'">4.0</xsl:when>
                <xsl:when test="$EdgeType='samples_cluster_composition'">3.0</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='Color'>
            <xsl:choose>
                <xsl:when test="$EdgeType='bait_prey'">#BB8888<xsl:value-of select="$HexIntensity"/></xsl:when>
                <xsl:when test="$EdgeType='bimap_block'">#5588BB<xsl:value-of select="$HexIntensity"/></xsl:when>
                <xsl:when test="$EdgeType='samples_cluster_composition'">#993300</xsl:when>
            </xsl:choose>
        </xsl:variable>
        <xsl:variable name='SourceArrow' select="none"/>
        <xsl:variable name='TargetArrow'>
            <xsl:choose>
                <xsl:when test="$EdgeType='bait_prey'">standard</xsl:when>
                <xsl:when test="$EdgeType='bimap_block'">standard</xsl:when>
                <xsl:otherwise>none</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>
        <xsl:copy use-attribute-sets="id edge_nodes">
            <xsl:attribute name='directed'><xsl:value-of select='$IsDirected'/></xsl:attribute>
            <xsl:copy-of select='g:data'/>
            <data key="y:eg">
                <y:PolyLineEdge>
                    <y:LineStyle color="{$Color}" type="{$LineStyle}"
                        width="{$Width}"/>
                    <xsl:if test="boolean($IsDirected)">
                        <y:Arrows source="{$SourceArrow}" target="{$TargetArrow}"/>
                    </xsl:if>
                </y:PolyLineEdge>
            </data>
        </xsl:copy>
    </xsl:template>
</xsl:stylesheet>
