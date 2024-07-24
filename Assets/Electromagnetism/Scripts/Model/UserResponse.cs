using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public class UserResponse
{
    public string name;
    public List<string> responses = new List<string>();

    public UserResponse(List<string> responses, string name)
    {
        this.responses = responses;
        this.name = name;
    }
}
