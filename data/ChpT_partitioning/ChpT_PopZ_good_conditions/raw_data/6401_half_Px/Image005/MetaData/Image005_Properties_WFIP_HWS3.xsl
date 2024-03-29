<?xml version="1.0"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
]>
<xsl:stylesheet version="1.1" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:dt="urn:schemas-microsoft-com:datatypes">
	<xsl:output method="html"/>
	<xsl:param name="tempVal" select="none"/>
	<xsl:template match="/">
		<HTML>
			<HEAD>
				<TITLE>
					<xsl:value-of select="@Name"/>
				</TITLE>
				<script type="text/javascript">
					  var timestampVisible = 0;
					  var LocalizeGSDEventsDetailsVisible = 0;
					  var GISTDrawHighresImageDetailsVisible = 0;
					  var GSDMergeEventListVisible = 0;
					  var FilterGSDEventListVisible = 0;
					  var DriftCompensationVisible = 0;
					  
					  function Show()
					  {
						  if(timestampVisible == 0)
						  {
							  window.document.getElementById("ID_1").style.display = "none";
							  window.document.getElementById("ID_2").style.display = "block";
							  timestampVisible = 1;
						  }
						  else
						  {
							  window.document.getElementById("ID_1").style.display = "block";
							  window.document.getElementById("ID_2").style.display = "none";
							  timestampVisible = 0;
						  }
					  }	
					  
					  function ShowLocalizeGSDEventsDetails()
					  {
						  if(LocalizeGSDEventsDetailsVisible == 0)
						  {
							  var a = window.document.getElementsByName("ID_2001");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "none";
							  }
					  
							  var b = window.document.getElementsByName("ID_2002");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "block";
							  }							  
							  
							  LocalizeGSDEventsDetailsVisible = 1;
						  }
						  else
						  {
							  var a = window.document.getElementsByName("ID_2001");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "block";
							  }
							  
							  var b = window.document.getElementsByName("ID_2002");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "none";
							  }							  
							  
							  LocalizeGSDEventsDetailsVisible = 0;
						  }
					  }
					  
					  function ShowGISTDrawHighresImageDetails()
					  {
					   if(GISTDrawHighresImageDetailsVisible == 0)
						  {
							  var a = window.document.getElementsByName("ID_2003");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "none";
							  }
							  
							  var b = window.document.getElementsByName("ID_2004");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "block";
							  }							  
							  
							  GISTDrawHighresImageDetailsVisible = 1;
						  }
						  else
						  {
							  var a = window.document.getElementsByName("ID_2003");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "block";
							  }
							  
							  var b = window.document.getElementsByName("ID_2004");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "none";
							  }							  
							  
							  GISTDrawHighresImageDetailsVisible = 0;
						  }
					  }
					  
					  function ShowGSDMergeEventListDetails()
					  {
					   if(GSDMergeEventListVisible == 0)
						  {
							  var a = window.document.getElementsByName("ID_2005");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "none";
							  }
							  
							  var b = window.document.getElementsByName("ID_2006");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "block";
							  }							  
							  
							  GSDMergeEventListVisible = 1;
						  }
						  else
						  {
							  var a = window.document.getElementsByName("ID_2005");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "block";
							  }
							  
							  var b = window.document.getElementsByName("ID_2006");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "none";
							  }							  
							  
							  GSDMergeEventListVisible = 0;
						  }
					  }
					  
					  function ShowFilterGSDEventListDetails()
					  {
					  	  if(FilterGSDEventListVisible == 0)
						  {
							  var a = window.document.getElementsByName("ID_2007");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "none";
							  }
							  
							  var b = window.document.getElementsByName("ID_2008");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "block";
							  }							  
							  
							  FilterGSDEventListVisible = 1;
						  }
						  else
						  {
							  var a = window.document.getElementsByName("ID_2007");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "block";
							  }
							  
							  var b = window.document.getElementsByName("ID_2008");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "none";
							  }							  
							  
							  FilterGSDEventListVisible = 0;
						  }
					  }
					  
					  function ShowDriftCompensationDetails()
					  {
					  	  if(DriftCompensationVisible == 0)
						  {
							  var a = window.document.getElementsByName("ID_2009");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "none";
							  }
							  
							  var b = window.document.getElementsByName("ID_2010");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "block";
							  }							  
							  
							  DriftCompensationVisible = 1;
						  }
						  else
						  {
							  var a = window.document.getElementsByName("ID_2009");
							  for(var i = 0; i != a.length; i++)
							  {
								  a[i].style.display = "block";
							  }
							  
							  var b = window.document.getElementsByName("ID_2010");
							  for(var j = 0; j != b.length; j++)
							  {
								  b[j].style.display = "none";
							  }							  
							  
							  DriftCompensationVisible = 0;
						  }
					  }
			</script>
			</HEAD>
			<BODY topmargin="0px" leftmargin="0px" bgcolor="#EEEEEE">
				<xsl:apply-templates select="Data"/>
			</BODY>
		</HTML>
	</xsl:template>
	<xsl:template match="Data">
		<xsl:apply-templates select="Image/ImageDescription"/>
		<xsl:apply-templates select="Image/TimeStampList"/>
		<xsl:apply-templates select="Image/Attachment"/>
		<xsl:apply-templates select="GISTEventList/GISTEventListDescription"/>
	</xsl:template>
	<xsl:template match="ATLCameraSettingDefinition">
		<br/>
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<tr>
				<th>ATLCameraSettingDefinition </th>
			</tr>
			<TR>
				<TD>
					<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
							<TD>Attibute Name</TD>
							<TD>Attibutes Value</TD>
						</TR>
						<xsl:for-each select="@*">
							<xsl:sort select="name()"/>
							<TR>
								<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<xsl:value-of select="name()"/>
								</TD>
								<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<xsl:value-of select="."/>
								</TD>
							</TR>
						</xsl:for-each>
					</TABLE>
				</TD>
			</TR>
			<TR>
				<TD>
					<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<xsl:for-each select="child::*">
							<xsl:sort select="child::text()"/>
							<TR>
								<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
									<xsl:value-of select="name()"/>
								</TD>
								<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
									<xsl:value-of select="."/>
								</TD>
							</TR>
							<xsl:for-each select="@*">
								<xsl:sort select="name()"/>
								<TR>
									<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<xsl:value-of select="name()"/>
									</TD>
									<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<xsl:value-of select="."/>
									</TD>
								</TR>
							</xsl:for-each>
						</xsl:for-each>
					</TABLE>
				</TD>
			</TR>
			<tr>
				<td>
					<xsl:apply-templates/>
				</td>
			</tr>
		</TABLE>
	</xsl:template>
	<xsl:template match="User-Comment">
  </xsl:template>
	<xsl:template name="break">
		<xsl:param name="text" select="//User-Comment"/>
		<xsl:comment>This inserts line breaks into the user description in place of line feeds</xsl:comment>
		<xsl:choose>
			<xsl:when test="contains($text, '&#xa;')">
				<xsl:value-of select="substring-before($text, '&#xa;')"/>
				<br/>
				<xsl:call-template name="break">
					<xsl:with-param name="text" select="substring-after($text,'&#xa;')"/>
				</xsl:call-template>
			</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="$text"/>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:template>
	<xsl:template match="ImageDescription">
		<xsl:variable name="isFileSaved">
			<xsl:if test="FileLocation != ' '">1</xsl:if>
		</xsl:variable>
		<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
					<TABLE width="100%" align="center" border="0" cellspacing="0" cellpadding="3" bgcolor="#FFFFFF">
						<TR>
							<TD align="center" valign="center" bgcolor="#EEEEEE" >
								<xsl:choose>
									<xsl:when test="0 != string-length(@ThumbnailPNGImage)">
											<img style="align:center; valign:middle" border="0" alt="thumbnail" src="data:image/png;base64,{@ThumbnailPNGImage}" />
									</xsl:when>
									<xsl:otherwise>
									</xsl:otherwise>
								</xsl:choose>
							</TD>
							<TD>
								<TABLE width="100%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD height="20" width="35%">Image: </TD>
										<TD>
											<B>
												<xsl:value-of select="Name"/>
											</B>
										</TD>
									</TR>
									<xsl:variable name="experimentType">
										<xsl:comment>Variable zum Filtern des !-Zeichens</xsl:comment>
										<xsl:choose>
											<xsl:when test="contains(//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/@ScanMode, '!')">
												<xsl:value-of select="substring-before(//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/@ScanMode, '!')"/>
											</xsl:when>
											<xsl:otherwise>
												<xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/@ScanMode"/>
											</xsl:otherwise>
										</xsl:choose>
									</xsl:variable>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
										<TD width="35%">Experiment type : </TD>
										<TD>
											<xsl:value-of select="$experimentType"/>
										</TD>
									</TR>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/@UserManagementUserName != 'UserManagementFeatureInactive'">
                      <xsl:variable name="UserManagementUserName">
                        <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/@UserManagementUserName"/>
                      </xsl:variable>
                      <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                        <TD width="40%">
                          User :
                        </TD>
                        <TD>
                          <xsl:value-of select="$UserManagementUserName"/>
                        </TD>
                      </TR>
                    </xsl:when>
                  </xsl:choose>                  
									<xsl:if test="$isFileSaved = '1'">
										<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<TD width="35%">File location : </TD>
											<TD>
												<xsl:value-of select="FileLocation"/>
											</TD>
										</TR>
									</xsl:if>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">Size : </TD>
										<TD>
											<xsl:value-of select="Size"/>
										</TD>
									</TR>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">Start Time : </TD>
										<TD>
											<xsl:value-of select="StartTime"/>
										</TD>
									</TR>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">End Time : </TD>
										<TD>
											<xsl:value-of select="EndTime"/>
										</TD>
									</TR>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">Total Exposures : </TD>
										<TD>
											<xsl:value-of select="FrameCount"/>
										</TD>
									</TR>
									<xsl:if test="//RelativeFocusCorrection != ' '">
										<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<TD width="35%">Relative Focus Correction : </TD>
											<TD>
												<xsl:value-of select="RelativeFocusCorrection"/>
											</TD>
										</TR>
									</xsl:if>
								</TABLE>
							</TD>
							<TD align="center" valign="top" rowspan="2">
								<A href="http://www.confocal-microscopy.com/" target="about:blank">
									<IMG src="LeicaLogo.jpg" border="0" alt="Leica Microsystems Heidelberg GmbH"/>
								</A>
							</TD>
						</TR>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<xsl:if test="//User-Comment != ' '">
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
								<TD colspan="2" width="35%">
									<xsl:call-template name="break"/>
								</TD>
							</TR>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<BR/>
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
					<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
							<TD>Dimension</TD>
							<TD>Logical Size</TD>
							<TD>Physical Length</TD>
							<TD>Physical Origin</TD>
							<TD>Voxel Size</TD>
						</TR>
						<xsl:for-each select="Dimensions/DimensionDescription">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD>
									<xsl:value-of select="@DimID"/>
								</TD>
								<TD>
									<xsl:value-of select="@NumberOfElements"/> &nbsp; <xsl:value-of select="@LogicalUnit"/>
								</TD>
								<TD>
									<xsl:value-of select="@Length"/> &nbsp;<xsl:value-of select="@Unit"/>
								</TD>
								<TD>
									<xsl:value-of select="@Origin"/> &nbsp;
									<xsl:choose>
										<xsl:when test="@UnitOrigin != ''">
											<xsl:value-of select="@UnitOrigin"/>
										</xsl:when>
										<xsl:otherwise>
											<xsl:value-of select="@Unit"/>
										</xsl:otherwise>
									</xsl:choose>
								</TD>
								<TD>
									<xsl:choose>
										<xsl:when test="@Voxel != ''">
											<xsl:value-of select="format-number(@Voxel, '0.00000')"/> &nbsp;<xsl:value-of select="@Unit"/>
										</xsl:when>
										<xsl:otherwise> --- </xsl:otherwise>
									</xsl:choose>
								</TD>
							</TR>
						</xsl:for-each>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<BR/>
		<xsl:variable name="isStereoTLFinetuningBase">
			<xsl:comment>Variable for display Stereo TLBase Finetuning-Settings</xsl:comment>
			<xsl:for-each select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator">
				<xsl:value-of select="@IsStereoTLFinetuningBase"/>
			</xsl:for-each>
		</xsl:variable>
		<xsl:variable name="typeOfStereoTLBase">
			<xsl:comment>Variable for type of Stereo TLBase</xsl:comment>
			<xsl:for-each select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator">
				<xsl:value-of select="@TypeOfStereoTLBase"/>
			</xsl:for-each>
		</xsl:variable>
		<xsl:variable name="isObliqueMode">
			<xsl:comment>Variable for display Oblique mode</xsl:comment>
			<xsl:choose>
				<xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ObliqueMode = '0'">0</xsl:when>
				<xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ObliqueMode = '1'">1</xsl:when>
				<xsl:otherwise>0</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="ChannelCount">
			<xsl:comment>Variable for channel count</xsl:comment>
      <xsl:value-of select="//Data/Image/ImageDescription/NumberOfChannels"/>
		</xsl:variable>    
    <xsl:variable name="CurrentChannelIndex">
			<xsl:comment>Variable for current channel index</xsl:comment>
      <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@CurrentChannelIndex"/>
		</xsl:variable>    
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
				<TD colspan="2">
					<TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="0" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
							<TD width="20%">
								<b>Channel Name </b>
							</TD>
							<TD width="20%">
								<b>Name </b>
							</TD>
							<TD width="15%">
								<b>Cube </b>
							</TD>
							<TD width="15%">
								<b>Contrast Method </b>
							</TD>
							<xsl:if test="isObliqueMode != '0'">
							  <TD width="7%">
								  <b>Oblique Value </b>
							  </TD>
							</xsl:if>
							<xsl:if test="$typeOfStereoTLBase = '0' or $typeOfStereoTLBase = '-1'">
								<TD width="7%">
									<b>Intensity </b>
								</TD>
							</xsl:if>
							<xsl:if test="$typeOfStereoTLBase = '0'">
								<TD width="7%">
									<b>CCIC intensity </b>
								</TD>
							</xsl:if>
							<TD width="8%">
								<b>Peak Emission </b>
							</TD>
							<TD width="9%">
								<b>Peak Excitation </b>
							</TD>
						</TR>
						<xsl:for-each select="//ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo[@Channel &gt;= 2000 and @Channel &lt; 3000]">
							<xsl:variable name="MyCounter" select="position()-1"/>
             	<xsl:variable name="VisibleChannel">
			          <xsl:comment>Variable for channel visibility</xsl:comment>
                <xsl:choose>
                   <xsl:when test="$ChannelCount = '1'">
                      <xsl:choose>
                        <xsl:when test="$CurrentChannelIndex != $MyCounter">0</xsl:when>   
                        <xsl:otherwise>1</xsl:otherwise>                    
                      </xsl:choose> 
                   </xsl:when>
                  <xsl:otherwise>1</xsl:otherwise>                    
                </xsl:choose> 
		          </xsl:variable> 
              
              <xsl:if test="$VisibleChannel = '1'">
              
                <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
								  <TD>
                    Channel <xsl:value-of select="position()"/>
								  </TD>
								  <TD>
									  <xsl:value-of select="@UserDefName"/>
								  </TD>
								  <TD>
									  <xsl:value-of select="@FluoCubeName"/> &nbsp;
                  </TD>
								  <TD>
                    <xsl:choose>
                      <xsl:when test="@ContrastingMethodName = 'GIST'">dSTORM / GSDIM</xsl:when>
                      <xsl:otherwise>
									      <xsl:value-of select="@ContrastingMethodName"/>
                      </xsl:otherwise>
                    </xsl:choose>
								  </TD>
								  <xsl:if test="$typeOfStereoTLBase = '-1'">
									  <TD>
										  <xsl:value-of select="@Intensity"/>
									  </TD>
								  </xsl:if>
									<xsl:if test="$typeOfStereoTLBase = '0'">
										<TD>
											<xsl:value-of select="@TL_Light-Intensity"/>
										</TD>
									</xsl:if>
									<xsl:if test="$typeOfStereoTLBase = '0'">
									  <TD>
										  <xsl:value-of select="@CCICCurrentPos"/>
									  </TD>
								  </xsl:if>
								  <xsl:if test="isObliqueMode != '0'">
									  <TD>
										  <xsl:value-of select="@ObliqueValue"/>
									  </TD>
								  </xsl:if>
								  <TD>
									  <xsl:value-of select="//Data/Image/ImageDescription/FluoDescription/FluoDescriptionRecord[@Identifier = concat('Channel ', $MyCounter)]/@Emission"/> &nbsp; nm
                  </TD>
								  <TD>
                    <xsl:choose>
                    <xsl:when test="//Data/Image/ImageDescription/FluoDescription/FluoDescriptionRecord[@Identifier = concat('Channel ', $MyCounter)]/@Excitation > 0">
                      <xsl:value-of select="//Data/Image/ImageDescription/FluoDescription/FluoDescriptionRecord[@Identifier = concat('Channel ', $MyCounter)]/@Excitation"/> &nbsp; nm
                    </xsl:when>
                    <xsl:otherwise>
                      0 &nbsp; nm
                    </xsl:otherwise>
                  </xsl:choose>
                  </TD>
							  </TR>             
                
              </xsl:if>
						</xsl:for-each>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<BR/>
    
		<xsl:variable name="isColorCamera">
			<xsl:if test="//ATLCameraSettingDefinition/CameraFormat/@IsGray ='0'">1</xsl:if>
		</xsl:variable>
    
		<xsl:variable name="canDoEMGain">
			<xsl:comment>Variable for display  EM Gain Settings</xsl:comment>
			<xsl:if test="//ATLCameraSettingDefinition/@CanDoEMGain ='1'">1</xsl:if>
		</xsl:variable>
    
		<xsl:variable name="isStereoIris">
			<xsl:comment>Variable for display Stereo Iris-Settings</xsl:comment>
			<xsl:for-each select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator">
				<xsl:value-of select="@CanDoStereoIris"/>
			</xsl:for-each>
		</xsl:variable>

    <xsl:variable name="showDIC">
      <xsl:comment>Variable for display DIC</xsl:comment>
      <xsl:for-each select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator">
        <xsl:value-of select="@CanDoDIC"/>
      </xsl:for-each>
    </xsl:variable>

    <xsl:variable name="showCondenser">
      <xsl:comment>Variable for display Condenser</xsl:comment>
      <xsl:for-each select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator">
        <xsl:value-of select="@CanDoCondenser"/>
      </xsl:for-each>
    </xsl:variable>
    
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
				<TD colspan="2">
					<TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="0" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
							<TD width="15%">
								<b>Channel Name </b>
							</TD>
							<xsl:if test="$isColorCamera != '1'">
								<TD width="15%">
									<b>LUT Name </b>
								</TD>
							</xsl:if>
							<TD width="15%">
								<b>Exposure Time </b>
							</TD>
							<TD width="7%">
								<b>Gain </b>
							</TD>
							<xsl:if test="$canDoEMGain = '1'">
								<TD width="10%">
									<b>EM Gain</b>
								</TD>
							</xsl:if>
							<xsl:if test="$isColorCamera = '1'">
								<TD width="15%">
									<b>Color-Gain Mode</b>
								</TD>
							</xsl:if>
							<xsl:if test="$isStereoIris = '1'">
								<TD width="15%">
									<b>Iris</b>
								</TD>
							</xsl:if>
							<xsl:if test="$showDIC = '1'">
							  <TD width="7%">
								<b>DIC </b>
							  </TD>
							</xsl:if>
							<xsl:if test="$showCondenser = '1'">
							  <TD width="7%">
								<b>Condenser </b>
							  </TD>
							</xsl:if>
							<TD width="15%">
								<b>Resolution XY</b>
							</TD>
							<TD width="15%">
								<b>Resolution Z</b>
							</TD>
						</TR>						
					  <xsl:for-each select="//ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo[@Channel &gt;= 2000 and @Channel &lt; 3000]">
						<xsl:variable name="here">
						  <xsl:comment>Variable for channel index</xsl:comment>
						  <xsl:choose>
							<xsl:when test="$ChannelCount = '1'">
							  <xsl:value-of select="$CurrentChannelIndex + 1"/>
							</xsl:when>
							<xsl:otherwise>
							  <xsl:value-of select="position()"/>
							</xsl:otherwise>
						  </xsl:choose>
						</xsl:variable>

						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
						  <TD>
							Channel <xsl:value-of select="position()"/> <!--/ <xsl:value-of select="$ChannelCount"/> , <xsl:value-of select="$CurrentChannelIndex + 1"/>, <xsl:value-of select="$here"/>-->
						  </TD>
						  <xsl:if test="$isColorCamera != '1'">
							<TD>
							  <xsl:value-of select="//Data/Image/ImageDescription/Channels/ChannelDescription[position()=$here]/@LUTName"/>
							</TD>
						  </xsl:if>
						  <TD>
							<xsl:value-of select="@ExposureTime"/>
						  </TD>
						  <TD>
							<xsl:value-of select="@Gain"/>
						  </TD>
						  <xsl:if test="$canDoEMGain = '1'">
							<TD>
							  <xsl:value-of select="@EMGain"/>
							</TD>
						  </xsl:if>
						  <xsl:if test="$isColorCamera = '1'">
							<TD>
							  ( R:<xsl:value-of select="@ColorGainValueR"/> -
							  G:<xsl:value-of select="@ColorGainValueG"/> -
							  B:<xsl:value-of select="@ColorGainValueB"/> )
							</TD>
						  </xsl:if>
						  <xsl:if test="$isStereoIris = '1'">
							<TD>
							  <xsl:value-of select="@IrisValue"/>
							</TD>
						  </xsl:if>
						  <xsl:if test="$showDIC = '1'">
							<TD>
								<xsl:choose>
									<xsl:when test="@DICLensName != ''">
									  <xsl:value-of select="@DICLensName"/>
									</xsl:when>
									<xsl:otherwise> --- </xsl:otherwise>
								</xsl:choose>							  
							</TD>
						  </xsl:if>
						  <xsl:if test="$showCondenser = '1'">
							<TD>
							  <xsl:value-of select="@CondenserLensName"/>
							</TD>
						  </xsl:if>
						  <TD>
							<xsl:choose>
							  <xsl:when test="//Data/Image/ImageDescription/Channels/ChannelDescription[position()=$here]/@OpticalResolutionXY != ''">
								<xsl:value-of select="//Data/Image/ImageDescription/Channels/ChannelDescription[position()=$here]/@OpticalResolutionXY"/>
							  </xsl:when>
							  <xsl:otherwise> --- </xsl:otherwise>
							</xsl:choose>
						  </TD>
						  <TD>
							<xsl:choose>
							  <xsl:when test="//Data/Image/ImageDescription/Channels/ChannelDescription[position()=$here]/@OpticalResolutionZ != ''">
								<xsl:value-of select="//Data/Image/ImageDescription/Channels/ChannelDescription[position()=$here]/@OpticalResolutionZ"/>
							  </xsl:when>
							  <xsl:otherwise> --- </xsl:otherwise>
							</xsl:choose>
						  </TD>
						</TR>
					  </xsl:for-each>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<xsl:variable name="canDoFFW">
			<xsl:comment>Variable for display FFW Settings</xsl:comment>
			<xsl:for-each select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo">
				<xsl:if test="@FFW_Emission1FilterName != '' or 
                      @FFW_Emission2FilterName != '' or 
                      @FFW_Excitation1FilterName != '' or 
                      @FFW_Excitation2FilterName != ''">1        
        </xsl:if>
			</xsl:for-each>
		</xsl:variable>
		<xsl:if test="$canDoFFW != ''">
			<HR width="98%"/>
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>Fast Filter Wheels</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
					<TD colspan="2">
						<TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="0" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
								<TD width="20%">
									<b>Channel Name </b>
								</TD>
								<TD width="20%">
									<b>Emission1 </b>
								</TD>
								<TD width="20%">
									<b>Emission2</b>
								</TD>
								<TD width="20%">
									<b>Excitation1</b>
								</TD>
								<TD width="20%">
									<b>Excitation2 </b>
								</TD>
							</TR>
							<xsl:for-each select="//ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo[@Channel &gt;= 2000 and @Channel &lt; 3000]">
             	   <xsl:variable name="MyCounter" select="position()-1"/>
                 <xsl:variable name="VisibleChannel">
			            <xsl:comment>Variable for channel visibility</xsl:comment>
                  <xsl:choose>
                     <xsl:when test="$ChannelCount = '1'">
                        <xsl:choose>
                          <xsl:when test="$CurrentChannelIndex != $MyCounter">0</xsl:when>   
                          <xsl:otherwise>1</xsl:otherwise>                    
                        </xsl:choose> 
                     </xsl:when>
                    <xsl:otherwise>1</xsl:otherwise>                    
                  </xsl:choose> 
		            </xsl:variable>

                <xsl:if test="$VisibleChannel = '1'">

  							    <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
									    <TD>
                        Channel <xsl:value-of select="position()"/>
									    </TD>
									    <TD>
											<xsl:choose>
												<xsl:when test="@FFW_Emission1FilterName != ''">
													<xsl:value-of select="@FFW_Emission1FilterName"/>
												</xsl:when>
												<xsl:otherwise> --- </xsl:otherwise>
											</xsl:choose>										    
									    </TD>
									    <TD>
											<xsl:choose>
												<xsl:when test="@FFW_Emission2FilterName != ''">
													<xsl:value-of select="@FFW_Emission2FilterName"/>
												</xsl:when>
												<xsl:otherwise> --- </xsl:otherwise>
											</xsl:choose>
									    </TD>
									    <TD>
										    <xsl:choose>
												<xsl:when test="@FFW_Excitation1FilterName != ''">
													<xsl:value-of select="@FFW_Excitation1FilterName"/>
												</xsl:when>
												<xsl:otherwise> --- </xsl:otherwise>
											</xsl:choose>
									    </TD>
									    <TD>
										    <xsl:choose>
												<xsl:when test="@FFW_Excitation2FilterName != ''">
													<xsl:value-of select="@FFW_Excitation2FilterName"/>
												</xsl:when>
												<xsl:otherwise> --- </xsl:otherwise>
											</xsl:choose>											
									    </TD>
								    </TR>              
                
                </xsl:if>                 
							</xsl:for-each>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<xsl:variable name="isILLEDAvailable">
			<xsl:comment>Variable for ILLED-Settings available</xsl:comment>
			<xsl:if test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '0'">1</xsl:if>
		</xsl:variable>
		<xsl:if test="$isILLEDAvailable = '1'">
      <xsl:variable name="canDoILLEDCombiLightMode">
        <xsl:comment>Variable for exsistence of an combilight slider</xsl:comment>
        <xsl:if test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDCanDoCombiLightSource = '1'">1</xsl:if>
      </xsl:variable>
			<HR width="98%"/>
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>Multicolor-Illumination</b>
					</TD>
				</TR>
			</TABLE>
      <xsl:choose>
        <xsl:when test="$canDoILLEDCombiLightMode = 1">
			    <TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				    <TR>
					    <TD>
						    <TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							    <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								    <TD width="40%">IL-LED Manual Slider Pos.</TD>
								    <xsl:choose>
									    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDCodedSliderPos = '1'">
										    <TD>LED</TD>
									    </xsl:when>
									    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDCodedSliderPos = '2'">
										    <TD>Mixed mode ( Lamp + LED )</TD>
									    </xsl:when>
									    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDCodedSliderPos = '3'">
										    <TD>Lamp</TD>
									    </xsl:when>
									    <xsl:otherwise>
										    <TD>undefined</TD>
									    </xsl:otherwise>
								    </xsl:choose>
							    </TR>
						    </TABLE>
					    </TD>
				    </TR>
			    </TABLE>
        </xsl:when>
      </xsl:choose>
			<xsl:variable name="isILLEDActive">
				<xsl:comment>Variable for ILLED active</xsl:comment>
				<xsl:choose>
					<xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDCodedSliderPos = '1'">1</xsl:when>
					<xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDCodedSliderPos = '2'">1</xsl:when>
				</xsl:choose>
			</xsl:variable>
      <xsl:variable name="canDoILLEDChangeWavelength">
				<xsl:comment>Variable if ILLED can change Wavelength</xsl:comment>
				<xsl:choose>
					<xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDCanDoChangeWavelength = '1'">1</xsl:when>
				</xsl:choose>
			</xsl:variable>
      <xsl:variable name="ILLEDUIMode" select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDUIMode"/>
			<xsl:if test="$isILLEDActive = '1'">
				<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
					<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
						<TD>
							<TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="0" bgcolor="#FFFFFF">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
									<TD width="25%">
										<b>Channel Name</b>
									</TD>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '0'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.1</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength0 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength0"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '1'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.2</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength1 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength1"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '2'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.3</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength2 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength2"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '3'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.4</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength3 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength3"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '4'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.5</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength4 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength4"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
                    <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '5'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.6</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength5 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength5"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '6'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.7</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength6 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength6"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
                  <xsl:choose>
                    <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '7'">
                      <TD>
                        <xsl:choose>
                          <xsl:when test="$canDoILLEDChangeWavelength = 1 or $ILLEDUIMode = 3">Illum.Chan.8</xsl:when>
                            <xsl:otherwise>
										          <xsl:choose>
											          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength7 > '0'">
												          <b>
													          <xsl:value-of select="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDWavelength7"/>nm
                                  </b>
											          </xsl:when>
											          <xsl:otherwise>-</xsl:otherwise>
										          </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
									    </TD>
                    </xsl:when>
                  </xsl:choose>
								</TR>
								<xsl:for-each select="//ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo[@Channel &gt;= 2000 and @Channel &lt; 3000]">
                   <xsl:variable name="MyCounter" select="position()-1"/>
                   <xsl:variable name="VisibleChannel">
			              <xsl:comment>Variable for channel visibility</xsl:comment>
                    <xsl:choose>
                       <xsl:when test="$ChannelCount = '1'">
                          <xsl:choose>
                            <xsl:when test="$CurrentChannelIndex != $MyCounter">0</xsl:when>   
                            <xsl:otherwise>1</xsl:otherwise>                    
                          </xsl:choose> 
                       </xsl:when>
                      <xsl:otherwise>1</xsl:otherwise>                    
                    </xsl:choose> 
		              </xsl:variable>

                  <xsl:if test="$VisibleChannel = '1'">

                      <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
										    <TD>
                          Channel <xsl:value-of select="position()"/>
										    </TD>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '0'">
                            <TD>
											        <xsl:choose>
												        <xsl:when test="@ILLEDActiveState0 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState0 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity0"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity0"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength0"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength0"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '1'">
                            <TD>
											        <xsl:choose>
												        <xsl:when test="@ILLEDActiveState1 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState1 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity1"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity1"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength1"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength1"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '2'">
                            <TD>
											        <xsl:choose>
												        <xsl:when test="@ILLEDActiveState2 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState2 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity2"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity2"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength2"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength2"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '3'">
                            <TD>
											       <xsl:choose>
												        <xsl:when test="@ILLEDActiveState3 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState3 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity3"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity3"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength3"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength3"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '4'">
                            <TD>
											        <xsl:choose>
												        <xsl:when test="@ILLEDActiveState4 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState4 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity4"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity4"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength4"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength4"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '5'">
                            <TD>
											        <xsl:choose>
												        <xsl:when test="@ILLEDActiveState5 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState5 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity5"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity5"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength5"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength5"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '6'">
                            <TD>
											        <xsl:choose>
												        <xsl:when test="@ILLEDActiveState6 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState6 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity6"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity6"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength6"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength6"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
										    <xsl:choose>
                          <xsl:when test="//Data/Image/Attachment[@Name='HardwareSetting']/ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ILLEDMaxWavelengths > '7'">
                            <TD>
											        <xsl:choose>
												        <xsl:when test="@ILLEDActiveState7 = '0'">-</xsl:when>
												        <xsl:when test="@ILLEDActiveState7 > '0'">
                                  <xsl:choose>
                                    <xsl:when test="$ILLEDUIMode = 1">
                                      <xsl:value-of select="@ILLEDIntensity7"/>%
                                    </xsl:when>
                                    <xsl:when test="$ILLEDUIMode = 2 or $ILLEDUIMode = 3">
                                      <xsl:value-of select="@ILLEDDiscreteIntensity7"/>%
                                    </xsl:when>
                                  </xsl:choose>
                                </xsl:when>
												        <xsl:otherwise>-</xsl:otherwise>
											        </xsl:choose>
                              <xsl:choose>
                                <xsl:when test="$canDoILLEDChangeWavelength = 1">
                                  (<xsl:value-of select="@ILLEDWavelength7"/> nm)
                                </xsl:when>
                                <xsl:when test="$ILLEDUIMode = 3">
                                  (<xsl:value-of select="@ILLEDDiscreteWavelength7"/>nm)
                                </xsl:when>
                              </xsl:choose>
										        </TD>
                          </xsl:when>
                        </xsl:choose>
									    </TR>
                  </xsl:if>
								</xsl:for-each>
							</TABLE>
						</TD>
					</TR>
				</TABLE>
			</xsl:if>
		</xsl:if>
	</xsl:template>
	<xsl:template match="TimeStampList">
		<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
			<TR>
				<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
					<b>Time Stamps:</b> &nbsp;
        </TD>
			</TR>
		</TABLE>
		<DIV ID="ID_1" style="display:block;">
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<TD>
                  Frame &nbsp; (<a href="javascript:Show()">Show All</a>)
                </TD>
								<TD>Relative Time (s)</TD>
								<TD>Absolute Time (h:m:s.ms)</TD>
								<TD>Date</TD>
							</TR>
							<xsl:for-each select="TimeStamp">
								<xsl:if test="not(position()!=1 and position()!=last())">
									<TR style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD>
											<xsl:number value="position()" format="1 "/>
										</TD>
										<TD>
											<xsl:choose>
												<xsl:when test="@RelativeTime != ''">
													<xsl:value-of select="@RelativeTime"/>
												</xsl:when>
												<xsl:otherwise> --- </xsl:otherwise>
											</xsl:choose>
										</TD>
										<TD>
											<xsl:value-of select="@Time"/>.<xsl:value-of select="@MiliSeconds"/>
										</TD>
										<TD>
											<xsl:choose>
												<xsl:when test="@Date != ''">
													<xsl:value-of select="@Date"/>
												</xsl:when>
												<xsl:otherwise> --- </xsl:otherwise>
											</xsl:choose>
										</TD>
									</TR>
								</xsl:if>
							</xsl:for-each>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
			<BR/>
		</DIV>
		<DIV ID="ID_2" style="display:none;">
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<TD>
                  Frame &nbsp; (<a href="javascript:Show()">Show first + last</a>)
                </TD>
								<TD>Relative Time</TD>
								<TD>Absolute Time</TD>
								<TD>Date</TD>
							</TR>
							<xsl:for-each select="TimeStamp">
								<TR style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD>
										<xsl:number value="position()" format="1 "/>
									</TD>
									<TD>
										<xsl:value-of select="@RelativeTime"/>
									</TD>
									<TD>
										<xsl:value-of select="@Time"/>.<xsl:value-of select="@MiliSeconds"/>
									</TD>
									<TD>
										<xsl:value-of select="@Date"/>
									</TD>
								</TR>
							</xsl:for-each>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</DIV>
		<BR/>
	</xsl:template>
	<xsl:template match="Attachment[@Name='HardwareSetting']">
		<HR width="98%"/>
		<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
			<TR>
				<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
					<b>Camera Settings</b>
				</TD>
			</TR>
		</TABLE>
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
					<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                Camera
              </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@CameraName"/>
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
								Format
						  </TD>
						  <TD>
							  <xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/@BinningText != ''">
										<xsl:value-of select="ATLCameraSettingDefinition/@BinningText"/>
									</xsl:when>
									<xsl:otherwise> --- </xsl:otherwise>
								</xsl:choose>
						  </TD>
						</TR>
						<xsl:variable name="isMicroScanningAvailanle">
							<xsl:comment>Variable for Micro Scanning available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoMicroScanning = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isMicroScanningAvailanle = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Micro Scanning</TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/@MicroScanning"/>
								</TD>
							</TR>
						</xsl:if>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                Digitization
              </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/CameraFormat/@Resolution"/>  &nbsp; bits
              </TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                Gamma
              </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@Gamma"/>
							</TD>
						</TR>
					  <xsl:variable name="canDoSpeedAndPort">
						<xsl:comment>Variable for Quality Mode available</xsl:comment>
						<xsl:if test="ATLCameraSettingDefinition/@CanDoChangeSpeedAndPort = '1'">1</xsl:if>
					  </xsl:variable>
					  <xsl:if test="$canDoSpeedAndPort = '1'">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                Quality Mode
              </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@SpeedandPort"/>
							</TD>
						</TR>
					  </xsl:if>
					  <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<TD width="40%">
						  Image Flip X Axis
						</TD>
						<TD>
						  <xsl:if test="ATLCameraSettingDefinition/@MirrorX = '1'">Yes</xsl:if>
						  <xsl:if test="ATLCameraSettingDefinition/@MirrorX = '0'">No</xsl:if>						  
						</TD>
					  </TR>
					  <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<TD width="40%">
						  Image Flip Y Axis
						</TD>
						<TD>
						  <xsl:if test="ATLCameraSettingDefinition/@MirrorY = '1'">Yes</xsl:if>
						  <xsl:if test="ATLCameraSettingDefinition/@MirrorY = '0'">No</xsl:if>
						</TD>
					  </TR>
						<xsl:variable name="isDualLightModeAvailable">
							<xsl:comment>Variable for Dual Light Mode available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoDualLightMode = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isDualLightModeAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Dual light mode</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/@DualLightMode = '0'">
										<TD>NIR mode</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/@DualLightMode = '1'">
										<TD>Standard mode</TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Undefined</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
						</xsl:if>
            <xsl:variable name="isColorCaptureModeAvailable">
              <xsl:comment>Variable for Color Capture Mode available</xsl:comment>
              <xsl:if test="ATLCameraSettingDefinition/@CanDoColorCapture = '1'">1</xsl:if>
            </xsl:variable>
            <xsl:if test="$isColorCaptureModeAvailable = '1'">
						  <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							  <TD width="40%">Color Capture Mode</TD>
							  <xsl:choose>
								  <xsl:when test="ATLCameraSettingDefinition/CameraFormat/@ColorCaptureMode = '0'">
									  <TD>Composite</TD>
								  </xsl:when>
								  <xsl:when test="ATLCameraSettingDefinition/CameraFormat/@ColorCaptureMode = '1'">
									  <TD>R-Channel only</TD>
								  </xsl:when>
								  <xsl:when test="ATLCameraSettingDefinition/CameraFormat/@ColorCaptureMode = '2'">
									  <TD>G-Channel only</TD>
								  </xsl:when>
								  <xsl:when test="ATLCameraSettingDefinition/CameraFormat/@ColorCaptureMode = '3'">
									  <TD>B-Channel only</TD>
								  </xsl:when>
								  <xsl:otherwise>
									  <TD>Undefined</TD>
								  </xsl:otherwise>
							  </xsl:choose>
						  </TR>
            </xsl:if>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                Black-Value
              </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@BlackValue"/>
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                White-Value
              </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@WhiteValue"/>
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">Shading Correction</TD>
							<xsl:choose>
								<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@OnlineShadingCorrection = '0'">
									<TD>OFF</TD>
								</xsl:when>
								<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@OnlineShadingCorrection = '1'">
									<TD>ON</TD>
								</xsl:when>
								<xsl:otherwise>
									<TD>Undefined</TD>
								</xsl:otherwise>
							</xsl:choose>
						</TR>
						<xsl:for-each select="//Data/Image/Attachment[@Name='Lost Frame Info']/Info">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                  Lost Frame Detected
                </TD>
								<TD>
									<xsl:choose>
										<xsl:when test="@LostFrameDetected = '0'">NO</xsl:when>
										<xsl:when test="@LostFrameDetected = '1'">YES</xsl:when>
									</xsl:choose>
								</TD>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                  Lost Sequence Detected
                </TD>
								<TD>
									<xsl:choose>
										<xsl:when test="@LostSequenceDetected = '0'">NO</xsl:when>
										<xsl:when test="@LostSequenceDetected = '1'">YES</xsl:when>
									</xsl:choose>
								</TD>
							</TR>
						</xsl:for-each>
						<xsl:variable name="isLiveBinningAvailable">
							<xsl:comment>Variable for Live Format available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoLiveBinning = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isLiveBinningAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Live Format</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/@IsLiveBinningActive = '0'">
										<TD>OFF</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/@IsLiveBinningActive = '1'">
									  <TD>
										<xsl:value-of select="ATLCameraSettingDefinition/@LiveBinningText"/>
									  </TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Undefined</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
						</xsl:if>
						<xsl:variable name="isCoolerAvailable">
							<xsl:comment>Variable for Cooler available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoControlCooler = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isCoolerAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Cooler</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/@IsCoolerActive = '0'">
										<TD>OFF</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/@IsCoolerActive = '1'">
										<TD>ON</TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Undefined</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
						</xsl:if>
						<xsl:variable name="isFanAvailable">
							<xsl:comment>Variable for Fan available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoFanControl = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isFanAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Fan Control</TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/@FanControl"/>
								</TD>
							</TR>
						</xsl:if>
						<xsl:variable name="isControlTemperatureAvailable">
							<xsl:comment>Variable for Control Temperature available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoControlTemperature = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isControlTemperatureAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Temperature Control</TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/@TargetTemperature"/> Degrees [C°]
								</TD>
							</TR>
						</xsl:if>
						<xsl:variable name="isHotSpotCorrectionAvailable">
							<xsl:comment>Variable for Hot Spot Correction available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoHotSpotCorrection = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isHotSpotCorrectionAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Hotspot Correction Threshold</TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/@ThresholdHotspot"/>%</TD>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Hotspot Correction Exposure Start</TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/@ExposureStartTimeHotspot"/>
								</TD>
							</TR>
						</xsl:if>
						<xsl:variable name="isBrightnessCorrectionAvailable">
							<xsl:comment>Variable for Brightness Correction available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoBrightnessCorrection = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isBrightnessCorrectionAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Brightness Correction</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/@IsBrightnessCorrectionActive = '0'">
										<TD>OFF</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/@IsBrightnessCorrectionActive = '1'">
										<TD>ON</TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Undefined</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
						</xsl:if>
						<xsl:variable name="isColorSatAvailable">
							<xsl:comment>Variable for ColorSaturation available</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@CanDoColorSat = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$isColorSatAvailable = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Color Saturation</TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/@ColorSat"/>
								</TD>
							</TR>
						</xsl:if>
						<xsl:variable name="canDoEnhancedAcquisitionMode">
						  <xsl:comment>Variable for Enhanced Acqiosition Mode</xsl:comment>
						  <xsl:if test="ATLCameraSettingDefinition/@CanDoEnhancedAcquisitionMode = '1'">1</xsl:if>
						</xsl:variable>
						<xsl:if test="$canDoEnhancedAcquisitionMode = '1'">
						  <xsl:variable name="enhancedAcquisitionMode">
							<xsl:comment>Variable for Enhanced Acqiosition Mode</xsl:comment>
							<xsl:if test="ATLCameraSettingDefinition/@HDR_DigitalFusion_Mode = '0'">None</xsl:if>
							<xsl:if test="ATLCameraSettingDefinition/@HDR_DigitalFusion_Mode = '1'">HDR</xsl:if>
							<xsl:if test="ATLCameraSettingDefinition/@HDR_DigitalFusion_Mode = '2'">Digital Fusion</xsl:if>
						  </xsl:variable>
						  <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">Acquisition Mode</TD>
							<TD>
							  <xsl:value-of select="$enhancedAcquisitionMode"/>
							</TD>
						  </TR>
						</xsl:if>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<xsl:variable name="isACFActive">
			<xsl:comment>Variable for display Advanced-Camera-Features-Settings</xsl:comment>
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/AdvancedCameraFeatures/@Advanced_Camera_Features_Activated = '0'">0</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/AdvancedCameraFeatures/@Advanced_Camera_Features_Activated = '1'">1</xsl:when>
				<xsl:otherwise>0</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<xsl:if test="$isACFActive = '1'">
			<HR width="98%"/>
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>Advanced Camera Features</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
								<TD colspan="2">
                  Global Features
                </TD>
							</TR>
							<xsl:for-each select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/AdvancedCameraFeatures/CameraFeature">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">
										<xsl:value-of select="@Name"/>
									</TD>
									<xsl:variable name="ACFEnabled">
										<xsl:value-of select="@*[ contains(name(),'Enabled') ] "/>
									</xsl:variable>
									<xsl:choose>
										<xsl:when test="$ACFEnabled = '1'">
											<xsl:variable name="ACFValue">
												<xsl:value-of select="@*[ contains(name(),'Value') ] "/>
											</xsl:variable>
											<xsl:choose>
												<xsl:when test="$ACFValue = ''">
													<TD>
                                 Enabled
                               </TD>
												</xsl:when>
												<xsl:otherwise>
													<TD>
														<xsl:value-of select="$ACFValue"/>
													</TD>
												</xsl:otherwise>
											</xsl:choose>
										</xsl:when>
										<xsl:otherwise>
											<TD>
                             Disabled
                           </TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:for-each>
							<xsl:for-each select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/AdvancedCameraFeatures/*">
								<xsl:if test="contains(name(),'Ch')">
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
										<TD colspan="2">
                        Channel <xsl:value-of select="name()"/> Features
                      </TD>
									</TR>
									<xsl:for-each select="CameraFeature">
										<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<TD width="40%">
												<xsl:value-of select="@Name"/>
											</TD>
											<xsl:variable name="ACFEnabled">
												<xsl:value-of select="@*[ contains(name(),'Enabled') ] "/>
											</xsl:variable>
											<xsl:choose>
												<xsl:when test="$ACFEnabled = '1'">
													<xsl:variable name="ACFValue">
														<xsl:value-of select="@*[ contains(name(),'Value') ] "/>
													</xsl:variable>
													<xsl:choose>
														<xsl:when test="$ACFValue = ''">
															<TD>
                                Enabled
                              </TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>
																<xsl:value-of select="$ACFValue"/>
															</TD>
														</xsl:otherwise>
													</xsl:choose>
												</xsl:when>
												<xsl:otherwise>
													<TD>
                            Disabled
                          </TD>
												</xsl:otherwise>
											</xsl:choose>
										</TR>
									</xsl:for-each>
								</xsl:if>
							</xsl:for-each>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		
		<xsl:variable name="isStereoTLFinetuningBase">
			<xsl:comment>Variable for display Stereo TLBase Finetuning-Settings</xsl:comment>
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@IsStereoTLFinetuningBase = '0'">0</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@IsStereoTLFinetuningBase = '1'">1</xsl:when>
				<xsl:otherwise>0</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="typeOfStereoTLBase">
			<xsl:comment>Variable for type of Stereo TLBase</xsl:comment>
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TypeOfStereoTLBase = '0'">0</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TypeOfStereoTLBase = '1'">1</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TypeOfStereoTLBase = '2'">2</xsl:when>
				<xsl:otherwise>-1</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
   <xsl:variable name="canDoFREEIllumination">
      <xsl:comment>Variable for display Stereo TLBase Finetuning-Settings</xsl:comment>
      <xsl:choose>
        <xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@CanDoFREEIllumination = '0'">0</xsl:when>
        <xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@CanDoFREEIllumination = '1'">1</xsl:when>
        <xsl:otherwise>0</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="ShutterOrPowerState">
      <xsl:comment>In LAS X 1.6 and 1.7 only Industry applications showed "Power State On/Off instead of "Shutter State open/closed". Since 3.0 we always use Power State</xsl:comment>
        Power State
    </xsl:variable>
    <xsl:variable name="OpenOrOn">
        On
    </xsl:variable>
    <xsl:variable name="ClosedOrOff">
        Off
    </xsl:variable>

    <xsl:if test="($typeOfStereoTLBase = '1' or $typeOfStereoTLBase = '2') or $canDoFREEIllumination = '1'">
			<HR width="98%"/>
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>Illumination Settings</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<xsl:for-each select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo">
								<xsl:variable name="MyCounter" select="position()"/>
								<xsl:variable name="WFCPrefix">WFC</xsl:variable>
								<xsl:variable name="ChannelPrefix">
									<xsl:value-of select="concat($WFCPrefix, $MyCounter)"/>
								</xsl:variable>
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
									<TD>
                    Channel <xsl:value-of select="$MyCounter"/>
									</TD>
								</TR>
								<TR>
									<TD>
										<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
											<xsl:if test="@StereoTLZoomApertureControlMode and $typeOfStereoTLBase = '2'">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              TL Zoom Aperture Control Mode
                            </TD>
													<xsl:choose>
														<xsl:when test="@StereoTLZoomApertureControlMode = '0'">
															<TD>Off</TD>
														</xsl:when>
														<xsl:when test="@StereoTLZoomApertureControlMode = '1'">
															<TD>On</TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>Undefined</TD>
														</xsl:otherwise>
													</xsl:choose>
												</TR>
											</xsl:if>
											<xsl:if test="@StereoTLContrastMode and $typeOfStereoTLBase = '2'">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              TL Contrast Mode
                            </TD>
													<xsl:choose>
														<xsl:when test="@StereoTLContrastMode = '1'">
															<TD>TL-BF</TD>
														</xsl:when>
														<xsl:when test="@StereoTLContrastMode = '2'">
															<TD>TL-RC</TD>
														</xsl:when>
														<xsl:when test="@StereoTLContrastMode = '3'">
															<TD>TL-DF</TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>Undefined</TD>
														</xsl:otherwise>
													</xsl:choose>
												</TR>
											</xsl:if>
											<xsl:if test="@StereoTLBase_Position and ($typeOfStereoTLBase = '1' or $typeOfStereoTLBase = '2')">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              TL Position
                            </TD>
													<TD>
														<xsl:value-of select="@StereoTLBase_Position"/>
													</TD>
												</TR>
											</xsl:if>
											<xsl:if test="@StereoTLBase_Aperture and $typeOfStereoTLBase = '2'">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              TL Aperture
                            </TD>
													<TD>
														<xsl:value-of select="@StereoTLBase_Aperture"/>
													</TD>
												</TR>
											</xsl:if>
											<xsl:if test="@TL_Shutter">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              TL <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                            </TD>
													<xsl:choose>
														<xsl:when test="@TL_Shutter = '0'">
															<TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:when test="@TL_Shutter = '1'">
															<TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>Undefined</TD>
														</xsl:otherwise>
													</xsl:choose>
												</TR>
											</xsl:if>
											<xsl:if test="@TL_Light-Intensity">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              TL Light Intensity
                            </TD>
													<TD>
														<xsl:value-of select="@TL_Light-Intensity"/>
													</TD>
												</TR>
											</xsl:if>
											<xsl:if test="@IL_Shutter">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              IL <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                            </TD>
													<xsl:choose>
														<xsl:when test="@IL_Shutter = '0'">
															<TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:when test="@IL_Shutter = '1'">
															<TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>Undefined</TD>
														</xsl:otherwise>
													</xsl:choose>
												</TR>
											</xsl:if>
											<xsl:if test="@IL_Light-Intensity">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              IL Light Intensity
                            </TD>
													<TD>
														<xsl:value-of select="@IL_Light-Intensity"/>
													</TD>
												</TR>
											</xsl:if>
											<xsl:if test="@MCI_Shutter">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              MCI <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                            </TD>
													<xsl:choose>
														<xsl:when test="@MCI_Shutter = '0'">
															<TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:when test="@MCI_Shutter = '1'">
															<TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>Undefined</TD>
														</xsl:otherwise>
													</xsl:choose>
												</TR>
											</xsl:if>
											<xsl:if test="@MCI_Light-Intensity">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              MCI Light Intensity
                            </TD>
													<TD>
														<xsl:value-of select="@MCI_Light-Intensity"/>
													</TD>
												</TR>
											</xsl:if>
											<xsl:if test="@MCI_Light-Scene">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              MCI Light Scene
                            </TD>
													<TD>
														<xsl:value-of select="@MCI_Light-Scene"/>
													</TD>
												</TR>
											</xsl:if>
                      <!-- The ringlight shutter attribute has a different name for different ringlight light sources
                           For example RL_MEB110_Shutter for a MEB110 ringlight, or RL_DVM6_Shutter for the Sirius ringlight 
                           So we select the attributes which start with RL_ and contain _Shutter and take the first (and most likely only) one -->
											<xsl:if test="@*[starts-with(name(),'RL_') and contains(name(),'_Shutter')][1]">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              RL <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                            </TD>
													<xsl:choose>
														<xsl:when test="@*[starts-with(name(),'RL') and contains(name(),'_Shutter')][1] = '0'">
															<TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:when test="@*[starts-with(name(),'RL') and contains(name(),'_Shutter')][1] = '1'">
															<TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>Undefined</TD>
														</xsl:otherwise>
													</xsl:choose>
												</TR>
											</xsl:if>
                      <!-- The ringlight intensity attribute has a different name for different ringlight light sources
                           For example RL_MEB110_Light-Intensity for a MEB110 ringlight, or RL_DVM6_Light-Intensity for the Sirius ringlight 
                           So we select the attributes which start with RL_ and contain _Light-Intensity and take the first (and most likely only) one -->
											<xsl:if test="@*[starts-with(name(),'RL') and contains(name(),'_Light-Intensity')][1]">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              RL Intensity
                            </TD>
													<TD>
														<xsl:value-of select="@*[starts-with(name(),'RL') and contains(name(),'_Light-Intensity')][1]"/>
													</TD>
												</TR>
											</xsl:if>
                      <xsl:if test="@RL_Light-Adapter_Text">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            RL Light Adapter
                          </TD>
                          <TD>
                          	<xsl:value-of select="@RL_Light-Adapter"/>&nbsp;(<xsl:value-of select="@RL_Light-Adapter_Text"/>)
                          </TD>
                        </TR>
                      </xsl:if>
											<xsl:if test="@RL_Light-Scene">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              RL Light Scene
                            </TD>
													<TD>
														<xsl:value-of select="@RL_Light-Scene"/>
													</TD>
												</TR>
											</xsl:if>
											<xsl:if test="@SLI_Shutter">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              SLI <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                            </TD>
													<xsl:choose>
														<xsl:when test="@SLI_Shutter = '0'">
															<TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:when test="@SLI_Shutter = '1'">
															<TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
														</xsl:when>
														<xsl:otherwise>
															<TD>Undefined</TD>
														</xsl:otherwise>
													</xsl:choose>
												</TR>
											</xsl:if>
											<xsl:if test="@SLI_Light-Intensity">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              SLI Intensity
                            </TD>
													<TD>
														<xsl:value-of select="@SLI_Light-Intensity"/>
													</TD>
												</TR>
											</xsl:if>
											<xsl:if test="@SLI_Light-Scene">
												<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
													<TD width="40%">
                              SLI Light Scene
                            </TD>
													<TD>
														<xsl:value-of select="@SLI_Light-Scene"/>
													</TD>
												</TR>
											</xsl:if>
                      <xsl:if test="@HDI_Shutter">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            HDI <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                          </TD>
                          <xsl:choose>
                            <xsl:when test="@HDI_Shutter = '0'">
                              <TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:when test="@HDI_Shutter = '1'">
                              <TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:otherwise>
                              <TD>Undefined</TD>
                            </xsl:otherwise>
                          </xsl:choose>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@HDI_Light-Intensity">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            HDI Light Intensity
                          </TD>
                          <TD>
                            <xsl:value-of select="@HDI_Light-Intensity"/>
                          </TD>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@HDI_Light-Scene">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            HDI Light Scene
                          </TD>
                          <TD>
                            <xsl:value-of select="@HDI_Light-Scene"/>
                          </TD>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@CXI_Shutter">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            CXI <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                          </TD>
                          <xsl:choose>
                            <xsl:when test="@CXI_Shutter = '0'">
                              <TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:when test="@CXI_Shutter = '1'">
                              <TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:otherwise>
                              <TD>Undefined</TD>
                            </xsl:otherwise>
                          </xsl:choose>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@CXI_Light-Intensity">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            CXI Light Intensity
                          </TD>
                          <TD>
                            <xsl:value-of select="@CXI_Light-Intensity"/>
                          </TD>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@CXI_Light-Adapter_Text">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            CXI Light Adapter
                          </TD>
                          <TD>
														<xsl:value-of select="@CXI_Light-Adapter"/>&nbsp;(<xsl:value-of select="@CXI_Light-Adapter_Text"/>)
                          </TD>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@NVI_Shutter">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            NVI <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                          </TD>
                          <xsl:choose>
                            <xsl:when test="@NVI_Shutter = '0'">
                              <TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:when test="@NVI_Shutter = '1'">
                              <TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:otherwise>
                              <TD>Undefined</TD>
                            </xsl:otherwise>
                          </xsl:choose>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@NVI_Light-Intensity">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            NVI Light Intensity
                          </TD>
                          <TD>
                            <xsl:value-of select="@NVI_Light-Intensity"/>
                          </TD>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@BLI_Shutter">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            BLI <xsl:value-of select="$ShutterOrPowerState"></xsl:value-of>
                          </TD>
                          <xsl:choose>
                            <xsl:when test="@BLI_Shutter = '0'">
                              <TD>
                                <xsl:value-of select="$ClosedOrOff"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:when test="@BLI_Shutter = '1'">
                              <TD>
                                <xsl:value-of select="$OpenOrOn"></xsl:value-of>
                              </TD>
                            </xsl:when>
                            <xsl:otherwise>
                              <TD>Undefined</TD>
                            </xsl:otherwise>
                          </xsl:choose>
                        </TR>
                      </xsl:if>
                      <xsl:if test="@BLI_Light-Intensity">
                        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                          <TD width="40%">
                            BLI Light Intensity
                          </TD>
                          <TD>
                            <xsl:value-of select="@BLI_Light-Intensity"/>
                          </TD>
                        </TR>
                      </xsl:if>
                                          
										</TABLE>
									</TD>
								</TR>
							</xsl:for-each>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<xsl:variable name="isAutoFocusHSActive">
			<xsl:comment>Variable for display AutoFocusHS-Settings</xsl:comment>
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AutoFocusHSActive = '0'">0</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AutoFocusHSActive = '1'">1</xsl:when>
				<xsl:otherwise>0</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="ShowOnlyForAutoFocusHS">
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '0'">1</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '3'">1</xsl:when>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="ShowXYTSettingsForAutoFocusHS">
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '1'">1</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '3'">1</xsl:when>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="ShowAFCPositionForAutoFocusHSMode">
			<xsl:comment>Variable for display AFCPosition</xsl:comment>
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '1'">
					<xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSAFCOffset >= '0'">1</xsl:if>
				</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '3'">
					<xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSAFCOffset >= '0'">1</xsl:if>
				</xsl:when>
			</xsl:choose>
		</xsl:variable>
		<xsl:if test="$isAutoFocusHSActive = '1'">
			<HR width="98%"/>
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>Autofocus Settings</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Autofocus-Subsystem</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '0'">
										<TD>Highspeed Autofocus</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '1'">
										<TD>Adaptive Focus Control</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '2'">
										<TD>Continuous Adaptive Focus Control</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSSubsystem = '3'">
										<TD>Adaptive Focus Control + Highspeed Autofocus</TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Undefined</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
							<xsl:if test="$ShowAFCPositionForAutoFocusHSMode = '1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">AFC Position</TD>
									<TD>
										<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSAFCOffset"/>
									</TD>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowOnlyForAutoFocusHS='0'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Configuration-Mode</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSConfigMode = '0'">
											<TD>Userdefined</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSConfigMode = '1'">
											<TD>Optimized</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>Undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowOnlyForAutoFocusHS='1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Capture range</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSOperationMode = '0'">
											<TD>Local</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSOperationMode = '1'">
											<TD>Global</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Local range</TD>
									<TD>
										<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSRangeLocal"/>
                      microns
                    </TD>
								</TR>
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Global range</TD>
									<TD>
										<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSRangeGlobal"/>
                      microns
                    </TD>
								</TR>
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Precision</TD>
									<TD>
										<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSPrecision"/>
									</TD>
								</TR>
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Focus-Device</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSZUseMode = '1'">
											<TD>Galvo</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSZUseMode = '2'">
											<TD>Z-Wide</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSZUseMode = '4'">
											<TD>FineFocus</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowXYTSettingsForAutoFocusHS='1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Workflow - Timelapse</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowTimelapse = '1'">
											<TD>Execute first time only</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowTimelapse = '2'">
											<TD>Execute on all times</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowTimelapse = '3'">
											<TD>
												<xsl:text>Execute every </xsl:text>
												<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowTimelapseIterator"/>
												<xsl:text> times</xsl:text>
											</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowXYTSettingsForAutoFocusHS='1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Workflow - Stage positions</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowXY = '1'">
											<TD>Execute on first position only</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowXY = '2'">
											<TD>Execute on all positions</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowXY = '3'">
											<TD>
												<xsl:text>Execute every </xsl:text>
												<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSWorkflowXYIterator"/>
												<xsl:text> positions</xsl:text>
											</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowOnlyForAutoFocusHS='1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Use ROI</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSCameranRoiIsUsed = '0'">
											<TD>NO</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@AFHSCameranRoiIsUsed = '1'">
											<TD>YES</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<p> </p>
		<xsl:variable name="isBestFocusActive">
			<xsl:comment>Variable for display BestFocus-Settings</xsl:comment>
			<xsl:if test="ATLCameraSettingDefinition/@IsAutofocusOnStart='1'">1</xsl:if>
		</xsl:variable>
		<xsl:variable name="ShowOnlyForBestFocus">
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '0'">1</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '3'">1</xsl:when>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="ShowXYTSettingsForBestFocus">
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '1'">1</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '3'">1</xsl:when>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="ShowAFCPositionForBestFocusMode">
			<xsl:comment>Variable for display AFCPosition</xsl:comment>
			<xsl:choose>
				<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '1'">
					<xsl:if test="ATLCameraSettingDefinition/Autofocus-config/@AFCOffset >= '0'">1</xsl:if>
				</xsl:when>
				<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '3'">
					<xsl:if test="ATLCameraSettingDefinition/Autofocus-config/@AFCOffset >= '0'">1</xsl:if>
				</xsl:when>
			</xsl:choose>
		</xsl:variable>
		<xsl:if test="$isBestFocusActive = '1'">
			<HR width="98%"/>
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>Autofocus Settings</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">Autofocus-Subsystem</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '0'">
										<TD>Best Focus</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '1'">
										<TD>Adaptive Focus Control</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '2'">
										<TD>Continuous Adaptive Focus Control</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFSubsystem = '3'">
										<TD>Adaptive Focus Control + Best Focus</TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Undefined</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
							<xsl:value-of select="$ShowAFCPositionForBestFocusMode"/>
							<xsl:if test="$ShowAFCPositionForBestFocusMode = '1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">AFC Position</TD>
									<TD>
										<xsl:value-of select="ATLCameraSettingDefinition/Autofocus-config/@AFCOffset"/>
									</TD>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowOnlyForBestFocus = '1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Precision</TD>
									<TD>
										<xsl:value-of select="ATLCameraSettingDefinition/Autofocus-config/@Precision"/>
									</TD>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowOnlyForBestFocus = '1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Focus Range</TD>
									<TD>
										<xsl:value-of select="ATLCameraSettingDefinition/Autofocus-config/@FocusRange * 1000000"/>
                      microns
                    </TD>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowOnlyForBestFocus = '1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Focus-Device</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@ZUseMode = '1'">
											<TD>Galvo</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@ZUseMode = '2'">
											<TD>Z-Wide</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@ZUseMode = '4'">
											<TD>FineFocus</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowOnlyForBestFocus = '1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Analyze Type</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFAnalyseType = '1'">
											<TD>Confocal Optimized</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFAnalyseType = '2'">
											<TD>Intensity Based Method</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFAnalyseType = '3'">
											<TD>Noise Based Method</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFAnalyseType = '4'">
											<TD>Widefield Optimized</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@AFAnalyseType = '5'">
											<TD>Reflection Based Method</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowXYTSettingsForBestFocus='1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Workflow - Timelapse</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@WorkflowTimelapse = '1'">
											<TD>Execute first time only</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@WorkflowTimelapse = '2'">
											<TD>Execute on all times</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@WorkflowTimelapse = '3'">
											<TD>
												<xsl:text>Execute every </xsl:text>
												<xsl:value-of select="ATLCameraSettingDefinition/Autofocus-config/@WorkflowTimelapseIterator"/>
												<xsl:text> times</xsl:text>
											</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="$ShowXYTSettingsForBestFocus='1'">
								<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
									<TD width="40%">Workflow - Stage positions</TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@WorkflowXY = '1'">
											<TD>Execute on first position only</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@WorkflowXY = '2'">
											<TD>Execute on all positions</TD>
										</xsl:when>
										<xsl:when test="ATLCameraSettingDefinition/Autofocus-config/@WorkflowXY = '3'">
											<TD>
												<xsl:text>Execute every </xsl:text>
												<xsl:value-of select="ATLCameraSettingDefinition/Autofocus-config/@WorkflowXYIterator"/>
												<xsl:text> positions</xsl:text>
											</TD>
										</xsl:when>
										<xsl:otherwise>
											<TD>undefined</TD>
										</xsl:otherwise>
									</xsl:choose>
								</TR>
							</xsl:if>
							<xsl:if test="ATLCameraSettingDefinition/Autofocus-config/@UseFixSliceNumber='1'">
								<xsl:if test="$ShowOnlyForBestFocus = '1'">
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="40%">No. of Steps</TD>
										<TD>
											<xsl:value-of select="ATLCameraSettingDefinition/Autofocus-config/@FixSliceNumber"/>
										</TD>
									</TR>
								</xsl:if>
							</xsl:if>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<p> </p>
		<xsl:variable name="isTirf">
			<xsl:comment>Variable for display TIRF-Settings</xsl:comment>
			<xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@CanDoTIRF='1'">1</xsl:if>
		</xsl:variable>
		<xsl:variable name="isTirfChannelAvailable">
			<xsl:for-each select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo">
				<xsl:if test="@ContrastingMethod = '206' or @ContrastingMethod = '208'">
          1
        </xsl:if>
			</xsl:for-each>
		</xsl:variable>
		<xsl:variable name="TirfManualMode">
			<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TIRFManualMode"/>
		</xsl:variable>
		<xsl:variable name="ShowTIRFFastIlluminationSwitchingMode">
			<xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@IsTIRFFastIlluminationSwitchingActive='1'">1</xsl:if>
		</xsl:variable>
		<xsl:if test="$isTirf = '1' and $isTirfChannelAvailable != '' ">
			<HR width="98%"/>
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>TIRF Settings</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">TIRF Mode</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TIRFManualMode = '1'">
										<TD>Expert  Mode</TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Automatic Mode</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">TIRF Penetration Depth Mode</TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TIRFPenetrationDepthMode = '0'">
										<TD>Fast</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TIRFPenetrationDepthMode = '1'">
										<TD>Variable</TD>
									</xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@TIRFPenetrationDepthMode = '2'">
										<TD>Const</TD>
									</xsl:when>
									<xsl:otherwise>
										<TD>Undefined</TD>
									</xsl:otherwise>
								</xsl:choose>
							</TR>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="0" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
								<TD width="20%">
									<b>Channel Name</b>
								</TD>
								<TD>
									<b>Wavelength</b>
								</TD>
								<TD>
									<b>Laser Intensity</b>
								</TD>
								<TD>
									<b>PenetrationDepth</b>
								</TD>
								<TD>
									<xsl:choose>
										<xsl:when test="$TirfManualMode = '1'">
											<b>Aperture</b>
										</xsl:when>
										<xsl:otherwise>
											<b>PenetrationDepth-Index</b>
										</xsl:otherwise>
									</xsl:choose>
								</TD>
								<TD>
									<b>Azimuth</b>
								</TD>
								<xsl:if test="$ShowTIRFFastIlluminationSwitchingMode = '1'">
									<TD>
										<b>Fast Illum.Switching Mode</b>
									</TD>
								</xsl:if>
							</TR>
							<xsl:for-each select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/WideFieldChannelInfo">
								<xsl:if test="@Channel &gt;= '2000' and @Channel &lt; '3000'">
									<xsl:variable name="MyCounter" select="position()-1"/>
									<xsl:variable name="WFCPrefix">nWFC</xsl:variable>
									<xsl:variable name="WFCPostfix">TIRF</xsl:variable>
									<xsl:variable name="TIRFChannelPrefix">
										<xsl:value-of select="concat($WFCPrefix, $MyCounter, $WFCPostfix)"/>
									</xsl:variable>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
										<TD>
											<xsl:value-of select="concat('Channel ', @Channel - 1999)"/>
										</TD>
										<TD>
											<xsl:value-of select="@TIRF_WaveLength"/>nm
										</TD>
										<TD>
											<xsl:value-of select="@TIRF_LaserIntensity"/>%
										</TD>
										<TD>
											<xsl:value-of select="@TIRF_PenetrationDepth"/>nm
										</TD>
										<TD>
											<xsl:choose>
												<xsl:when test="$TirfManualMode = '1'">
													<xsl:value-of select="@TIRF_Aperture"/>
												</xsl:when>
												<xsl:otherwise>
													<xsl:value-of select="@TIRF_PenetrationDepthIndex"/>
												</xsl:otherwise>
											</xsl:choose>
										</TD>
										<TD>
											<xsl:value-of select="@TIRF_Azimuth"/>
										</TD>
										<xsl:if test="$ShowTIRFFastIlluminationSwitchingMode = '1'">
											<TD>
												<xsl:choose>
													<xsl:when test="@TIRF_FastIlluminationSwitchingMode = '1'">
                            FLUO
                          </xsl:when>
													<xsl:otherwise>
                            TIRF
                          </xsl:otherwise>
												</xsl:choose>
											</TD>
										</xsl:if>
									</TR>
								</xsl:if>
							</xsl:for-each>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<p> </p>
    <xsl:variable name="ShowStageRotationAngle">
      <xsl:comment>Only show the stage rotation angle if microscope is capable of stage rotation</xsl:comment>
      <xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@CanDoStageRotation != '0'">1</xsl:if>
    </xsl:variable>
    <xsl:variable name="ShowColumnTiltAngle">
      <xsl:comment>Only show the column tilt angle if microscope is capable of column tilting</xsl:comment>
      <xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@CanDoColumnTilt != '0'">1</xsl:if>
    </xsl:variable>
    <TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
			<TR>
				<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
					<b>Microscope Settings</b>
				</TD>
			</TR>
		</TABLE>
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
					<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                  System Name
                </TD>
							<TD>
								LAS X
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
								Microscope Model
							</TD>
							<TD>
									<xsl:choose>
										<xsl:when test="ATLCameraSettingDefinition/@MicroscopeModel != ''">
											<xsl:value-of select="ATLCameraSettingDefinition/@MicroscopeModel" />
										</xsl:when>
										<xsl:otherwise>
											<xsl:choose>
												<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@MicroscopeFamilyType = '0'">
												  Unknown
												</xsl:when>
												<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@MicroscopeFamilyType = '1'">
													<xsl:choose>
														<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@IsInverseMicroscopeModel != '1'">
															DMI4000-6000
														</xsl:when>
														<xsl:otherwise>
															DM5000-6000
														</xsl:otherwise>
													</xsl:choose>
												</xsl:when>
												<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@MicroscopeFamilyType = '2'">
													DMi8
												</xsl:when>
												<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@MicroscopeFamilyType = '3'">
													Stereo
												</xsl:when>
												<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@MicroscopeFamilyType = '4'">
													Manual
												</xsl:when>
												<xsl:otherwise>
													--
												</xsl:otherwise>
											</xsl:choose>
										</xsl:otherwise>
									</xsl:choose>
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                  Objective
                </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@ObjectiveName"/>
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                  Order number (Obj.)
                </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@ObjectiveNumber"/>
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                  Numerical aperture (Obj.)
                </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@NumericalAperture"/>
							</TD>
						</TR>
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                  Refraction index
                </TD>
							<TD>
								<xsl:value-of select="ATLCameraSettingDefinition/@RefractionIndex"/>
							</TD>
						</TR>
            <!-- Check whether there is a mag changer, only add the row, if there is one -->
            <xsl:choose>
              <xsl:when test="ATLCameraSettingDefinition/@TubeOpticName != ''">
                <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                  <TD width="40%">
                    Magnification-Changer
                  </TD>
                  <TD>
                    <xsl:value-of select="ATLCameraSettingDefinition/@TubeOpticName"/>x
                  </TD>
                </TR>
              </xsl:when>
            </xsl:choose>
            <!-- For Stereos: Show Zoom and Video magnification. Note that we only show the value of 'SameZoom', because the feature 'Same Zoom for each channel' is always on since LAS X 3.0 -->
            <xsl:choose>
              <xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@MicroscopeFamilyType = '3'">
                <xsl:variable name="ZoomValue">
                  <xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@SameZoom"/>
                </xsl:variable>
                <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                  <TD width="40%">
                    Zoom
                  </TD>
                  <TD>
                    <xsl:value-of select="format-number($ZoomValue,'0.##')"/>x
                  </TD>
                </TR>
                <xsl:choose>
                  <xsl:when test="ATLCameraSettingDefinition/@CanDoTotalVideoMag">
                    <xsl:variable name="TotalMagValue">
                      <xsl:value-of select="ATLCameraSettingDefinition/@TotalVideoMag"/>
                    </xsl:variable>
                    <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                      <TD width="40%">
                        Total Video Magnification
                      </TD>
                      <TD>
                        <xsl:value-of select="format-number($TotalMagValue,'0.#')"/>x
                      </TD>
                    </TR>
                  </xsl:when>
                </xsl:choose>
              </xsl:when>
            </xsl:choose>
            <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
							<TD width="40%">
                  Z Movement
                </TD>
							<TD>
								<xsl:choose>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@eStackFlowMode = '1'">
                      Lambda then Z
                    </xsl:when>
									<xsl:when test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@eStackFlowMode = '2'">
                      Z then Lambda
                    </xsl:when>
								</xsl:choose>
							</TD>
						</TR>
			<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
				<TD width="40%">
                  Camera Light
                </TD>
				<TD>
					<xsl:value-of select="ATLCameraSettingDefinition/@PortCameraPercentage"/>%
				</TD>
			</TR>
            <xsl:if test="$ShowColumnTiltAngle = '1'">
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Tilt Angle
                </TD>
                <TD>
                  <xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@ColumnTiltAngle"/>°
                </TD>
              </TR>
            </xsl:if>
						<xsl:if test="$ShowStageRotationAngle = '1'">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
	                  Stage Rotation Angle
	                </TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@StageRotationAngle"/>°
	              </TD>
							</TR>
						</xsl:if>
					</TABLE>
					<p> </p>
				</TD>
			</TR>
		</TABLE>
		<p> </p>
		<xsl:variable name="ShowDiaphragmSettings">
			<xsl:comment>Variable for display Diaphragm-Settings</xsl:comment>
			<xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@bIsDiaphragmStoreModeEnabled ='1'">1</xsl:if>
		</xsl:variable>
		<xsl:if test="$ShowDiaphragmSettings = '1'">
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>Diaphragm Settings</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                    IL-Field Diaphragm Settings
                  </TD>
								<TD>
									<xsl:value-of disable-output-escaping="yes" select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@FieldDiaphragmILSettingsToDisplay"/>
								</TD>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                    TL-Field Diaphragm Settings
                  </TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@m_sFieldDiaphragmTLSettings"/>
								</TD>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                    TL-Aperture Diaphragm Settings
                  </TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@m_sApertureDiaphragmTLSettings"/>
								</TD>
							</TR>
						</TABLE>
						<p> </p>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<p> </p>
		<xsl:variable name="ShowDICSettings">
			<xsl:comment>Variable for display DIC-Settings</xsl:comment>
			<xsl:if test="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@bIsDICStoreModeEnabled ='1'">1</xsl:if>
		</xsl:variable>
		<xsl:if test="$ShowDICSettings = '1'">
			<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
				<TR>
					<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
						<b>DIC Settings</b>
					</TD>
				</TR>
			</TABLE>
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                  DIC Prism
                </TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@PrismNameDIC"/>
								</TD>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                  Condenser Prism
                </TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@PrismNameCond"/>
								</TD>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD width="40%">
                  DIC Finetuning
                </TD>
								<TD>
									<xsl:value-of select="ATLCameraSettingDefinition/WideFieldChannelConfigurator/@DICBias"/>
								</TD>
							</TR>
						</TABLE>
						<p> </p>
					</TD>
				</TR>
			</TABLE>
    </xsl:if>

    
    <!-- The following block about FRAP WF is obsolete, because the FRAP WF informations are saved in a new created attachment since LAS X 2.0.
         I've not deleted it, because of the evaluation of experiments created with LAS X 1.5.0 in higher releases. -->
    <xsl:variable name="WFFrapExperiment">
      <xsl:value-of select="ATLCameraSettingDefinition/Block_FRAPWF/@FrapWFCanDoFrap"/>
    </xsl:variable>

    <xsl:if test="$WFFrapExperiment = '1'">

      <TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
        <TR>
          <TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
            <b>WF FRAP Settings</b>
          </TD>
        </TR>
      </TABLE>

      <TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
        <TR>
          <TD>
            <TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Diaphragm
                </TD>
                <TD>
                  <xsl:value-of select="ATLCameraSettingDefinition/Block_FRAPWF/@FrapWFDiaphragmPositionDescription"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser intensity
                </TD>
                <TD>
                  <xsl:value-of select="ATLCameraSettingDefinition/Block_FRAPWF/@FrapWFConvertedIntensityString"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser impulse length
                </TD>
                <TD>
                  <xsl:value-of select="ATLCameraSettingDefinition/Block_FRAPWF/@FrapWFLaserImpulseLength"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser impulse count
                </TD>
                <TD>
                  <xsl:value-of select="ATLCameraSettingDefinition/Block_FRAPWF/@FrapWFLaserImpulseNumber+1"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser impulse interim time
                </TD>
                <TD>
                  <xsl:value-of select="ATLCameraSettingDefinition/Block_FRAPWF/@FrapWFLaserImpulseInterimTime"/>
                </TD>
              </TR>
          </TABLE>
          </TD>
        </TR>
      </TABLE>
    </xsl:if>
      
    <xsl:variable name="climateControlsAvailable">
      <xsl:value-of select="ATLCameraSettingDefinition/ClimateControl/@ClimateControlCurrentlyInUse"/> 
    </xsl:variable>
    <xsl:if test="$climateControlsAvailable = '1'">

      <TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
        <TR>
          <TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
            <b>Environmental Settings</b>
          </TD>
        </TR>
      </TABLE>

      <TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
        <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
          <TD colspan="2">
            <TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="0" bgcolor="#FFFFFF">
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
                <TD width="20%">
                  <b>Unit Name </b>
                </TD>
                <TD width="10%">
                  <b>Enabled </b>
                </TD>
                <TD width="12%">
                  <b>Setpoint </b>
                </TD>
                <TD width="12%">
                  <b>Actual Value </b>
                </TD>
                <TD width="10%">
                  <b>Range Control </b>
                </TD>
                <TD width="12%">
                  <b>Minimum Range </b>
                </TD>
                <TD width="12%">
                  <b>Maximum Range </b>
                </TD>
                <TD width="12%">
                  <b>Stop On Alert </b>
                </TD>
              </TR>
              <xsl:for-each select="//ATLCameraSettingDefinition/ClimateControl/ClimateUnit">
                <xsl:variable name="MyCounter" select="position()-1"/>
                <xsl:variable name="SetpointVal" select="concat(@SetpointValue, '.0')"/>
                <xsl:variable name="SetpointValPointNull" select="@SetpointValue"/>

                <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
                  <TD>
                    <xsl:value-of select="@Name"/>
                  </TD>
                  <TD>
                    <xsl:if test="@Enabled != '0'">
                      ON
                    </xsl:if>
                    <xsl:if test="@Enabled = '0'">
                      OFF
                    </xsl:if>
                  </TD>
                  <TD>
                    <xsl:if test="@SetpointValue*10 mod 10 != 0">
                      <xsl:value-of select="@SetpointValue"/>
                    </xsl:if>
                    <xsl:if test="@SetpointValue*10 mod 10 = 0">
                      <xsl:value-of select="concat(@SetpointValue, '.0')"/>
                    </xsl:if>
                    <xsl:value-of select="@DimensionUnit"/>
                  </TD>
                  <TD>
                    <xsl:if test="@ActualValue*10 mod 10 != 0">
                      <xsl:value-of select="@ActualValue"/>
                    </xsl:if>
                    <xsl:if test="@ActualValue*10 mod 10 = 0">
                      <xsl:value-of select="concat(@ActualValue, '.0')"/>
                    </xsl:if>
                    <xsl:value-of select="@DimensionUnit"/>
                  </TD>
                  <TD>
                    <xsl:if test="@RangeControlEnabled != '0'">
                      ON
                    </xsl:if>
                    <xsl:if test="@RangeControlEnabled = '0'">
                      OFF
                    </xsl:if>
                  </TD>
                  <xsl:if test="@RangeControlEnabled != '0'">
                    <TD>
                      <xsl:if test="@RangeMinimumValue*10 mod 10 != 0">
                        <xsl:value-of select="@RangeMinimumValue"/>
                      </xsl:if>
                      <xsl:if test="@RangeMinimumValue*10 mod 10 = 0">
                        <xsl:value-of select="concat(@RangeMinimumValue, '.0')"/>
                      </xsl:if>
                      <xsl:value-of select="@DimensionUnit"/>
                    </TD>
                    <TD>
                      <xsl:if test="@RangeMaximumValue*10 mod 10 != 0">
                        <xsl:value-of select="@RangeMaximumValue"/>
                      </xsl:if>
                      <xsl:if test="@RangeMaximumValue*10 mod 10 = 0">
                        <xsl:value-of select="concat(@RangeMaximumValue, '.0')"/>
                      </xsl:if>
                      <xsl:value-of select="@DimensionUnit"/>
                    </TD>
                    <TD>
                      <xsl:if test="@RangeControlStopExperimentOnAlert != '0'">
                        ON
                      </xsl:if>
                      <xsl:if test="@RangeControlStopExperimentOnAlert = '0'">
                        OFF
                      </xsl:if>
                    </TD>
                  </xsl:if>
                  <xsl:if test="@RangeControlEnabled = '0'">
                    <TD>
                      -
                    </TD>
                    <TD>
                      -
                    </TD>
                    <TD>
                      -
                    </TD>
                  </xsl:if>
                </TR>
              </xsl:for-each>
            </TABLE>
          </TD>
        </TR>
      </TABLE>
      </xsl:if>
     
	</xsl:template>
	<xsl:template match="Attachment[@Name='ProcessingHistory']">
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5">
			<TR>
				<TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
					<B> Processing History: </B>
				</TD>
			</TR>
			<xsl:for-each select="ProcessingHistory">
				<xsl:sort select="@CreationTime" order="descending"/>
				<xsl:choose>
					<xsl:when test="@ProcessingType = 'OnlineDyeSeparation_Target' and contains(@TargetUniqueIDList, /Data/Image/ImageDescription/UniqueID ) ">
						<tr>
							<td>
								<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
									<TR>
										<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
											<B>
												<xsl:value-of select="@CreationTime"/>
											</B>
										</TD>
										<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
											<B> Applied OnlineDyeSeparation parameters </B>
										</TD>
									</TR>
								</TABLE>
								<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
									<TR>
										<TD>
											<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
												<TR style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
													<td>
														<xsl:copy-of select="TR"/>
													</td>
												</TR>
											</TABLE>
										</TD>
									</TR>
								</TABLE>
							</td>
						</tr>
						<tr>
							<td colspan="2" align="center" border="0" cellspacing="5" cellpadding="5"> &nbsp;</td>
						</tr>
					</xsl:when>
					<xsl:when test="@ProcessingType = 'OnlineDyeSeparation_Source' and contains(@SourceUniqueIDList, /Data/Image/ImageDescription/UniqueID ) ">
						<TR bgcolor="#DDDAD7">
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B>
									<xsl:value-of select="@CreationTime"/>
								</B>
							</TD>
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B> Active OnlineDyeSeparation parameters </B>
							</TD>
						</TR>
						<TR align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
							<TD>
								<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<TR style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>
											<xsl:copy-of select="TR"/>
										</td>
									</TR>
								</TABLE>
							</TD>
						</TR>
						<tr>
							<td colspan="2" align="center" border="0" cellspacing="5" cellpadding="5"> &nbsp;</td>
						</tr>
					</xsl:when>
					<xsl:when test="@ProcessingType = 'LocalizeGSDEvents' ">
						<TR border="1" bgcolor="#DDDAD7">
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B>
									<xsl:value-of select="@CreationTime"/>
								</B>
							</TD>
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B> GSD event localisation</B>
							</TD>
						</TR>
						<TR border="1" bgcolor="#DDDAD7">
							<TD colspan="2" border="0" cellspacing="0" cellpadding="0">
								<TABLE ID="ID_2001" style="display:block;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source Image</td>
										<td>
											<xsl:copy-of select="Data/Image/ImageDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>Source Image &nbsp;  (<a href="javascript:ShowLocalizeGSDEventsDetails()">Show Information</a>) </B>
										</td>
									</tr>
								</TABLE>
								<TABLE ID="ID_2002" style="display:none;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>Source Image &nbsp;  (<a href="javascript:ShowLocalizeGSDEventsDetails()">Hide</a>) </B>
										</td>
									</tr>
									<tr>
										<td colspan="2">
											<xsl:apply-templates select="Data"/>
										</td>
									</tr>
								</TABLE>
							</TD>
						</TR>
					</xsl:when>
					<xsl:when test="@ProcessingType = 'GISTDrawHighresImage' ">
						<TR topmargin="0" leftmargin="0" border="0" bgcolor="#DDDAD7">
							<TD topmargin="0" leftmargin="0" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<B>
									<xsl:value-of select="@CreationTime"/>
								</B>
							</TD>
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B> Draw GSD highresolution image</B>
							</TD>
						</TR>
						<TR border="1" bgcolor="#DDDAD7">
							<TD colspan="2" border="0" cellspacing="0" cellpadding="0">
								<TABLE ID="ID_2003" style="display:block;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>Source Image: &nbsp;  (<a href="javascript:ShowGISTDrawHighresImageDetails()">Show Information</a>) </B>
										</td>
									</tr>
								</TABLE>
								<TABLE ID="ID_2004" style="display:none;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>Source Image: &nbsp;  (<a href="javascript:ShowGISTDrawHighresImageDetails()">Hide</a>) </B>
										</td>
									</tr>
									<tr>
										<td colspan="2">
											<xsl:apply-templates select="Data"/>
										</td>
									</tr>
								</TABLE>
							</TD>
						</TR>
					</xsl:when>
					<xsl:when test="@ProcessingType = 'GSDMergeEventList' ">
						<TR topmargin="0" leftmargin="0" border="0" bgcolor="#DDDAD7">
							<TD topmargin="0" leftmargin="0" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<B>
									<xsl:value-of select="@CreationTime"/>
								</B>
							</TD>
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B> Merge GSD Eventlist</B>
							</TD>
						</TR>
						<TR border="1" bgcolor="#DDDAD7">
							<TD colspan="2" border="0" cellspacing="0" cellpadding="0">
								<TABLE ID="ID_2005" style="display:block;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>GSD Eventlist &nbsp;  (<a href="javascript:ShowGSDMergeEventListDetails()">Show Information</a>) </B>
										</td>
									</tr>
								</TABLE>
								<TABLE ID="ID_2006" style="display:none;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>GSD Eventlist &nbsp;  (<a href="javascript:ShowGSDMergeEventListDetails()">Hide</a>) </B>
										</td>
									</tr>
									<tr>
										<td colspan="2">
											<xsl:apply-templates select="Data"/>
										</td>
									</tr>
								</TABLE>
							</TD>
						</TR>
					</xsl:when>
					<xsl:when test="@ProcessingType = 'FilterGSDEventList' ">
						<TR topmargin="0" leftmargin="0" border="0" bgcolor="#DDDAD7">
							<TD topmargin="0" leftmargin="0" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<B>
									<xsl:value-of select="@CreationTime"/>
								</B>
							</TD>
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B> Filter GSD Eventlist</B>
							</TD>
						</TR>
						<TR border="1" bgcolor="#DDDAD7">
							<TD colspan="2" border="0" cellspacing="0" cellpadding="0">
								<TABLE ID="ID_2007" style="display:block;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>GSD Eventlist &nbsp;  (<a href="javascript:ShowFilterGSDEventListDetails()">Show Information</a>) </B>
										</td>
									</tr>
								</TABLE>
								<TABLE ID="ID_2008" style="display:none;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>GSD Eventlist &nbsp;  (<a href="javascript:ShowFilterGSDEventListDetails()">Hide</a>) </B>
										</td>
									</tr>
									<tr>
										<td colspan="2">
											<xsl:apply-templates select="Data"/>
										</td>
									</tr>
								</TABLE>
							</TD>
						</TR>
					</xsl:when>
					<xsl:when test="@ProcessingType = 'DriftCompensation' ">
						<TR topmargin="0" leftmargin="0" border="0" bgcolor="#DDDAD7">
							<TD topmargin="0" leftmargin="0" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<B>
									<xsl:value-of select="@CreationTime"/>
								</B>
							</TD>
							<TD align="left" style="font-family: arial, helvetica; font-size: 7pt; font-weight: bold; color: 000000; padding: 3px;">
								<B> Drift Compensation </B>
							</TD>
						</TR>
						<TR border="1" bgcolor="#DDDAD7">
							<TD colspan="2" border="0" cellspacing="0" cellpadding="0">
								<TABLE ID="ID_2009" style="display:block;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>GSD Eventlist &nbsp;  (<a href="javascript:ShowDriftCompensationDetails()">Show Information</a>) </B>
										</td>
									</tr>
								</TABLE>
								<TABLE ID="ID_2010" style="display:none;" topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<xsl:copy-of select="TR"/>
									<tr style="font-family: arial, helvetica; font-size: 7pt; font-weight: normal; color: 000000; padding: 3px;">
										<td>Name Of Source GSD EventList</td>
										<td>
											<xsl:copy-of select="Data/GISTEventList/GISTEventListDescription/Name"/>
										</td>                 
									</tr>
									<tr bgcolor="#DDDAD7">
										<td colspan="2" align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
											<B>GSD Eventlist &nbsp;  (<a href="javascript:ShowDriftCompensationDetails()">Hide</a>) </B>
										</td>
									</tr>
									<tr>
										<td colspan="2">
											<xsl:apply-templates select="Data"/>
										</td>
									</tr>
								</TABLE>
							</TD>
						</TR>
					</xsl:when>
				</xsl:choose>
			</xsl:for-each>
		</TABLE>
	</xsl:template>
	<xsl:template match="GISTEventListDescription">
		<TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
					<TABLE width="100%" align="center" border="0" cellspacing="0" cellpadding="3" bgcolor="#FFFFFF">
						<TR>
							<TD>
								<TABLE width="100%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD height="20" width="35%">EventList :</TD>
										<TD>
											<B>
												<xsl:value-of select="Name"/>
											</B>
										</TD>
									</TR>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">File name :</TD>
										<TD>
											<xsl:value-of select="FileLocation"/>
										</TD>
									</TR>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">Size :</TD>
										<TD>
											<xsl:value-of select="Size"/>
										</TD>
									</TR>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">CreationInfo :</TD>
										<TD>
											<xsl:for-each select="CreationInfo">
												<xsl:value-of select="@Date"/>
											</xsl:for-each>
										</TD>
									</TR>
									<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
										<TD width="35%">Number of Events :</TD>
										<TD style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
											<xsl:for-each select="NumberOfEvents">
												<xsl:value-of select="@NumberOfEventsValue"/>
											</xsl:for-each>
										</TD>
									</TR>
								</TABLE>
							</TD>
							<TD align="center" valign="top" rowspan="2">
								<A href="http://www.confocal-microscopy.com/" target="about:blank">
									<IMG src="LeicaLogo.jpg" border="0" alt="Leica Microsystems Heidelberg GmbH"/>
								</A>
							</TD>
						</TR>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<xsl:if test="//User-Comment != ' '">
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="center" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; padding: 3px;">
								<TD colspan="2" width="35%">
									<xsl:call-template name="break"/>
								</TD>
							</TR>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
		</xsl:if>
		<BR/>
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
					<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
							<TD>Threshold</TD>
							<TD>Gainfactor</TD>
							<TD>FieldOfView</TD>
						</TR>
						<xsl:for-each select="LocalizationParameters">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD>
									<xsl:value-of select="@Threshold"/>
								</TD>
								<TD>
									<xsl:value-of select="@Gain"/>
								</TD>
								<TD>
                  (<xsl:value-of select="@FieldOfViewX1"/> pixel, <xsl:value-of select="@FieldOfViewX1"/> pixel, <xsl:value-of select="@FieldOfViewX2"/> pixel, <xsl:value-of select="@FieldOfViewX2"/> pixel) (x1, y1, x2, y2)
              </TD>
							</TR>
						</xsl:for-each>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<BR/>
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
					<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
							<TD>SeriesID</TD>
							<TD>Source Image Name</TD>
							<TD>Experiment Path</TD>
						</TR>
						<xsl:for-each select="SourceData/GISTSourceDataDescription">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD>
									<xsl:value-of select="@SeriesID"/>
								</TD>
								<TD>
									<xsl:value-of select="@ImageName"/>
								</TD>
								<TD>
									<xsl:value-of select="@ExperimentPath"/>
								</TD>
							</TR>
						</xsl:for-each>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
		<BR/>
		<xsl:for-each select="DataAnalysis/XML3DCalibration/calibrationGSD3Ddata">
			<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
				<TR>
					<TD>
						<TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
								<TD>3D Calibration Name</TD>
								<TD>WaveLenght</TD>
							</TR>
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD>
									<xsl:value-of select="name"/>
								</TD>
								<TD>
									<xsl:value-of select="DyeWaveLengthInSample"/> nm
					</TD>
							</TR>
						</TABLE>
					</TD>
				</TR>
			</TABLE>
			<BR/>
		</xsl:for-each>
		<TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
			<TR>
				<TD>
          Eventlist export format:
          <TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
						<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: bold; color: 000000; padding: 3px;">
							<TD>Name</TD>
							<TD>Size in Bytes</TD>
							<TD>ColumnType</TD>
							<TD>Description</TD>
						</TR>
						<xsl:for-each select="Columns/GISTColumnsDescription">
							<TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
								<TD>
									<xsl:value-of select="@Name"/>
								</TD>
								<TD>
									<xsl:value-of select="@SizeInByte"/>
								</TD>
								<TD>
									<xsl:if test="@ColumnType = '1'">INT</xsl:if>
									<xsl:if test="@ColumnType = '2'">UINT</xsl:if>
									<xsl:if test="@ColumnType = '3'">FLOAT</xsl:if>
									<xsl:if test="@ColumnType = '4'">STRING</xsl:if>
								</TD>
								<TD>
									<xsl:if test="@Name = 'Dim4'">reserved for future use</xsl:if>
									<xsl:if test="@Name = 'frameID'">id of frame in t-series</xsl:if>
									<xsl:if test="@Name = 'eventD'">id for event inside a frame</xsl:if>
									<xsl:if test="@Name = 'channel1'">photon count of event in channel one</xsl:if>
									<xsl:if test="@Name = 'x1'">X position of event in channel one in pixel </xsl:if>
									<xsl:if test="@Name = 'y1'">Y position of event in channel one in pixel </xsl:if>
									<xsl:if test="@Name = 'z1'">Z position of event in channel one in meter </xsl:if>
									<xsl:if test="@Name = 'SeriesID'">the reference to the <EM>Source Image Name</EM> shown above</xsl:if>
									<xsl:if test="@Name = 'sigmaX'"> standard deviation in X direction </xsl:if>
									<xsl:if test="@Name = 'sigmaY'"> standard deviation in Y direction </xsl:if>
								</TD>
							</TR>
						</xsl:for-each>
					</TABLE>
				</TD>
			</TR>
		</TABLE>
	</xsl:template>
  <xsl:template match="Attachment[@Name='InteractiveMeasurement']"/>
  
  <xsl:template match="Attachment[@Name='ContextDescription']">
    
    <xsl:variable name="WFFrapExperiment">
      <xsl:value-of select="Wizard/Block_FRAPWF/@FrapWFCanDoFrap"/>
    </xsl:variable>

    <xsl:if test="$WFFrapExperiment = '1'">

      <TABLE width="98%" align="center" border="0" cellspacing="5" cellpadding="5">
        <TR>
          <TD align="left" style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
            <b>WF FRAP Settings</b>
          </TD>
        </TR>
      </TABLE>

      <TABLE width="98%" align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="#DDDAD7">
        <TR>
          <TD>
            <TABLE topmargin="0" leftmargin="0" width="100%" align="left" border="1" cellspacing="0" cellpadding="5" bgcolor="#FFFFFF">
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Diaphragm
                </TD>
                <TD>
                  <xsl:value-of select="Wizard/Block_FRAPWF/@FrapWFDiaphragmPositionDescription"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser intensity
                </TD>
                <TD>
                  <xsl:value-of select="Wizard/Block_FRAPWF/@FrapWFConvertedIntensityString"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser impulse length
                </TD>
                <TD>
                  <xsl:value-of select="Wizard/Block_FRAPWF/@FrapWFLaserImpulseLength"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser impulse count
                </TD>
                <TD>
                  <xsl:value-of select="Wizard/Block_FRAPWF/@FrapWFLaserImpulseNumber+1"/>
                </TD>
              </TR>
              <TR style="font-family: arial, helvetica; font-size: 8pt; font-weight: normal; color: 000000; padding: 3px;">
                <TD width="40%">
                  Laser impulse interim time
                </TD>
                <TD>
                  <xsl:value-of select="Wizard/Block_FRAPWF/@FrapWFLaserImpulseInterimTime"/>
                </TD>
              </TR>
          </TABLE>
          </TD>
        </TR>
      </TABLE>
    </xsl:if>
  
  </xsl:template>
  
  <xsl:template match="Attachment[@Name='Annotation']"/>

</xsl:stylesheet>
