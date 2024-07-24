// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

//It belongs to the Custom directory and is named GouraudShader
Shader "Custom/PhongShader"
{
	Properties
	{
		_MainColor("Color", Color) = (1, 1, 1, 1)
		//I defined the position of a light source myself
		_LightPos("LightPosition", Vector) = (-0.5, 5.5, -0.5, 1)
	}

		SubShader
	{
		Tags{ "RenderType" = "Opaque" }
		LOD 100

		Pass
	{
		CGPROGRAM
#pragma vertex vert
#pragma fragment frag

#include "UnityCG.cginc"

		struct v2f
	{
		float4 vertex : SV_POSITION;
		//TEXCOORD is used for interpolation
		float3 worldPos : TEXCOODR0;
		float3 worldNor : TEXCOODR1;
	};

	float4 _MainColor;
	float4 _LightPos;

	v2f vert(float4 vertex : POSITION, float3 normal : NORMAL)
	{
		v2f o;
		o.vertex = UnityObjectToClipPos(vertex);

		//Convert the coordinate information and normal information of the point we handwritten before to the world coordinate system for calculation.
		//Interpolate this information to the frag function
		o.worldPos = mul(UNITY_MATRIX_M, vertex);
		o.worldNor = mul(UNITY_MATRIX_M, normal);

		return o;
	}

	//The point of Phong Shading is to calculate in the frag function
	fixed4 frag(v2f i) : SV_Target
	{
		//Direction of light
		float3 lightDir = normalize(_LightPos - i.worldPos);
		//Distance, used to calculate the attenuation of light
		float3 dist = distance(_LightPos, i.worldPos);

		//Calculation of Diffuse Light.
		float lightPor = max(0, dot(i.worldNor, lightDir));
		//Attenuation coefficient, divided by 2 is a mysterious reason
		float atten = 2 / pow(dist, 2);

		//Return the final calculation result
		return _MainColor * lightPor * atten;
	}
		ENDCG
	}
	}
}