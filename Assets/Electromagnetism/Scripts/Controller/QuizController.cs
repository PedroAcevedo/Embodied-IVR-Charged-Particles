using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class QuizController : MonoBehaviour
{
    public GameObject instructions;
    public GameObject backToSimulation;
    public GameObject floatingFigure;
    public GameObject recordingUI;
    public GameObject optionPanel;
    public GameObject recordingIndicator;
    public GameObject previousButton;
    public GameObject finalUI;
    public List<GameObject> buttonList;
    public List<Button> MenuButtons;
    public TMPro.TextMeshProUGUI questionLabel;
    public TMPro.TextMeshProUGUI currentQuestion;
    public TMPro.TextMeshProUGUI recordButtonLabel;

    public static Questionnaire questionnaire;


    private string selectedOption = "";
    private List<string> questions;
    public List<int> isRandomized;
    private int questionCounter = 0;
    private int currentScene = 0;
    private bool isRecording = false;
    private MicrophoneController microphone;

    // Start is called before the first frame update
    void Start()
    {
        Debug.Log(Application.persistentDataPath);

        for (int i = 0; i < buttonList.Count; i++)
        {
            setButtonEvent(i);
        }

        questions = new List<string>();
        isRandomized = new List<int>();

        var jsonText = Resources.Load<TextAsset>("Data/experimentation");

        questionnaire = JsonUtility.FromJson<Questionnaire>(jsonText.text);
        updateCurrentQuestionLabel();
        setCurrentQuestion();
        buttonState();

        microphone = new MicrophoneController(GetComponent<AudioSource>());
    }

    // Update is called once per frame
    void Update()
    {

    }

    public void onOptionSelected(int position)
    {
        resetCurrentSelection();

        selectedOption = questionnaire.scenes[currentScene].questions[questionCounter].options[position];
        questionnaire.scenes[currentScene].questions[questionCounter].selectedAnswer = selectedOption;
        buttonList[position].GetComponent<UnityEngine.UI.Outline>().enabled = true;

    }

    public void changeQuestion(int direction)
    {
        if (isRecording) stopRecording();

        resetCurrentSelection();

        questionCounter = questionCounter + direction;

        if (questionCounter == (questionnaire.scenes[currentScene].questions.Count))
        {
            instructions.SetActive(false);
            floatingFigure.SetActive(false);
            backToSimulation.SetActive(true);

            questionCounter = questionCounter - 1;
        }
        else
        {
            updateCurrentQuestionLabel();
            setCurrentQuestion();
            buttonState();
        }
    }

    public void recordAnswer()
    {
        if (!isRecording)
        {
            isRecording = true;
            microphone.startRecording();
            recordButtonLabel.text = "Stop Recording";
            recordingIndicator.SetActive(true);
        }
        else
        {
            stopRecording();
        }

    }

    public void returnToSurvey()
    {
        backToSimulation.SetActive(false);
        floatingFigure.SetActive(true);
        instructions.SetActive(true);

        setSelectedAnswer();
    }

    public void changeScene()
    {
        currentScene = SimulationController.currentScene;
        questionCounter = 0;
        isRandomized.Clear();

        if (currentScene == 4)
        {
            finalUI.SetActive(true);
        }
        else
        {
            returnToSurvey();
            updateCurrentQuestionLabel();
            setCurrentQuestion();
            buttonState();
        }
    }

    private void setButtonEvent(int position)
    {
        buttonList[position].GetComponent<UnityEngine.UI.Button>().onClick.AddListener(delegate {
            onOptionSelected(position);
        });
    }

    private void resetCurrentSelection()
    {
        if (selectedOption != "" && questionnaire.scenes[currentScene].questions[questionCounter].options.Count > 0)
        {
            if (isAnswerInOptions() != -1)
            {
                Debug.Log("Current Selection -> " + selectedOption);
                Debug.Log("Current Selection Index -> " + isAnswerInOptions());

                buttonList[questionnaire.scenes[currentScene].questions[questionCounter].options.IndexOf(selectedOption)].GetComponent<UnityEngine.UI.Outline>().enabled = false;
            }
        }
    }

    private void updateCurrentQuestionLabel()
    {
        string question = "Question #" + (questionCounter + 1);

        currentQuestion.text = question;
    }

    private void setCurrentQuestion()
    {
        questionLabel.text = questionnaire.scenes[currentScene].questions[questionCounter].body;
        
        floatingFigure.SetActive(questionnaire.scenes[currentScene].questions[questionCounter].hasImage);
        previousButton.SetActive(questionCounter != 0);

        if (floatingFigure.activeSelf)
        {
            floatingFigure.transform.GetChild(0).GetChild(0).GetComponent<Image>().sprite = loadImageFromResources(questionnaire.scenes[currentScene].name, questionnaire.scenes[currentScene].questions[questionCounter].figure);
        }

        if (questionnaire.scenes[currentScene].questions[questionCounter].options.Count > 0)
        {
            optionPanel.SetActive(true);
            recordingUI.SetActive(false);

            if (!isRandomized.Contains(questionCounter))
            {
                isRandomized.Add(questionCounter);

                questionnaire.scenes[currentScene].questions[questionCounter].options = fisherYatesShuffle(questionnaire.scenes[currentScene].questions[questionCounter].options);
            }


            if (!questionnaire.scenes[currentScene].questions[questionCounter].hasImage)
            {
                loadOptionImages(questionnaire.scenes[currentScene].name);
            }
            else
            {
                loadLabelOptions();
            }

            setSelectedAnswer();
        }
        else
        {
            optionPanel.SetActive(false);
            recordingUI.SetActive(true);
        }
    }

    private Sprite loadImageFromResources(string folder, string image)
    {
        return Resources.Load<Sprite>("Questionnaire/" + folder + "/" + image);
    }

    private void loadOptionImages(string folder)
    {
        for(int i = 0; i < buttonList.Count; i++)
        {
            buttonList[i].transform.GetChild(0).gameObject.SetActive(false);
            buttonList[i].transform.GetChild(1).gameObject.SetActive(true);

            buttonList[i].transform.GetChild(1).GetComponent<Image>().sprite = loadImageFromResources(folder, questionnaire.scenes[currentScene].questions[questionCounter].options[i]);
        }
    }

    private void loadLabelOptions()
    {
        for (int i = 0; i < buttonList.Count; i++)
        {
            buttonList[i].transform.GetChild(0).gameObject.SetActive(true);
            buttonList[i].transform.GetChild(1).gameObject.SetActive(false);

            buttonList[i].transform.GetChild(0).GetComponent<TMPro.TextMeshProUGUI>().text = questionnaire.scenes[currentScene].questions[questionCounter].options[i];
        }
    }

    private void buttonState()
    {
        MenuButtons[0].interactable = true;
        MenuButtons[1].interactable = true;

        if (questionCounter == 0)
        {
            MenuButtons[0].interactable = false;
        }

        //if (questionCounter == (questionnarie.scenes[currentScene].questions.Count - 1))
        //{
        //    MenuButtons[1].interactable = false;
        //}
    }

    private void setSelectedAnswer()
    {
        Debug.Log("Set selected answer ->" + questionnaire.scenes[currentScene].questions[questionCounter].selectedAnswer);

        if (questionnaire.scenes[currentScene].questions[questionCounter].selectedAnswer != "")
        {
            selectedOption = questionnaire.scenes[currentScene].questions[questionCounter].selectedAnswer;
            buttonList[questionnaire.scenes[currentScene].questions[questionCounter].options.IndexOf(selectedOption)].GetComponent<UnityEngine.UI.Outline>().enabled = true;
        }
    }

    private int isAnswerInOptions()
    {
        return questionnaire.scenes[currentScene].questions[questionCounter].options.IndexOf(selectedOption);
    }

    private void stopRecording()
    {
        isRecording = false;
        recordButtonLabel.text = "Start Recording";
        recordingIndicator.SetActive(false);
        microphone.stopRecording(SimulationController.playerID + "_scene_" + (currentScene + 1) + "_question_" + (questionCounter + 1) + "_");
        questionnaire.scenes[currentScene].questions[questionCounter].recordingFile = microphone.filename;
    }

    private List<string> fisherYatesShuffle(List<string> list)
    {
        List<string> tempList = new List<string>(list);

        for (int i = 0; i < list.Count; i++)
        {
            int j = Random.Range(0, (list.Count - 1));

            string tempOpt = tempList[i];
            tempList[i] = tempList[j];
            tempList[j] = tempOpt;
        }

        return tempList;
    }
}
