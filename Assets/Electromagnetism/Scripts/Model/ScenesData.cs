using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public class ScenesData
{
    public string userId;
    public bool isRightHanded;
    public List<SceneData> scenes = new List<SceneData>();
    public List<UserResponse> experimentation = new List<UserResponse>();
    public string sceneRecording;

    public ScenesData(string userId, bool isRightHanded, string sceneRecording, List<SceneData> scenes, List<UserResponse> experimentation)
    {
        this.userId = userId;
        this.sceneRecording = sceneRecording;
        this.isRightHanded = isRightHanded;
        this.scenes = scenes;
        this.experimentation = experimentation;
    }
}
