#include "NaoWidget.hpp"

#include "dart/external/imgui/imgui.h"

#include "NaoWorldNode.hpp"

//==============================================================================
NaoWidget::NaoWidget(
    dart::gui::osg::ImGuiViewer* viewer,
    NaoWorldNode* node,
    dart::simulation::WorldPtr world)
  : mViewer(viewer),
    mNode(node),
    mWorld(world),
    mGuiVelGain(10),
    mGuiZmpGain(10),
    mGuiComHeight(0.33),
    mGuiReferenceVelocityX(0.0),
    mGuiReferenceVelocityY(0.0),
    mGuiReferenceVelocityOmega(0.0),
    mComHeight(mGuiComHeight),
    mGuiHeadlights(true),
    mGuiControlMode(2),
    mGuiBeheavior(0),
    mControlMode(2),
    mBeheavior(0),
    mExternalForceMode(0),
    mGuiExternalForceMode(0),
    mGuiExternalForceStartFrame(250),
    mGuiExternalForceX(0.0),
    mGuiExternalForceY(0.0),
    mGuiExternalForceZ(0.0),
    mGuiExternalForcePeriodicPhase(0.0),
    mGuiExternalForcePeriodicFrequency(1.0),
    mGuiComRoll(mNode->getController()->getBalanceFootPos()(0)),
    mGuiComPitch(mNode->getController()->getBalanceFootPos()(1)),
    mGuiComYaw(mNode->getController()->getBalanceFootPos()(2)),
    mGuiComX(mNode->getController()->getBalanceFootPos()(3)),
    mGuiComY(mNode->getController()->getBalanceFootPos()(4)),
    mGuiComZ(mNode->getController()->getBalanceFootPos()(5))
{
  // Do nothing
}

//==============================================================================
void NaoWidget::render()
{
  ImGui::SetNextWindowPos(ImVec2(10,20));
  if (!ImGui::Begin("Nao Control", nullptr, ImVec2(360,640), 0.5f,
                    ImGuiWindowFlags_MenuBar |
                    ImGuiWindowFlags_HorizontalScrollbar))
  {
    // Early out if the window is collapsed, as an optimization.
    ImGui::End();
    return;
  }

  // Menu
  if (ImGui::BeginMenuBar())
  {
    if (ImGui::BeginMenu("Menu"))
    {
      if (ImGui::MenuItem("Exit"))
        mViewer->setDone(true);
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Help"))
    {
      if (ImGui::MenuItem("About DART"))
        mViewer->showAbout();
      ImGui::EndMenu();
    }
    ImGui::EndMenuBar();
  }

  ImGui::Text("NAO Robot Controlled by MPC");
  ImGui::Spacing();

  if (ImGui::CollapsingHeader("Help"))
  {
    ImGui::PushTextWrapPos(ImGui::GetCursorPos().x + 320);
    ImGui::Text("User Guide:\n");
    ImGui::Text("%s", mViewer->getInstructions().c_str());
    ImGui::Text("Press [r] to reset Nao to the initial position.\n");
    ImGui::Text("Press [a] to push forward Nao torso.\n");
    ImGui::Text("Press [s] to push backward Nao torso.\n");
    ImGui::Text("Press [d] to push left Nao torso.\n");
    ImGui::Text("Press [f] to push right Nao torso.\n");
    ImGui::Text("Left-click on a block to select it.\n");
    ImGui::PopTextWrapPos();
  }

  if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen))
  {
    int e = mViewer->isSimulating() ? 0 : 1;
    if(mViewer->isAllowingSimulation())
    {
      if (ImGui::RadioButton("Play", &e, 0) && !mViewer->isSimulating())
        mViewer->simulate(true);
      ImGui::SameLine();
      if (ImGui::RadioButton("Pause", &e, 1) && mViewer->isSimulating())
        mViewer->simulate(false);
    }

    // Balance Point
    ImGui::RadioButton("Beheavior: Walk", &mGuiBeheavior, 1);
    ImGui::RadioButton("Beheavior: Balance", &mGuiBeheavior, 0);
    ImGui::RadioButton("Beheavior: Free Balance", &mGuiBeheavior, 2);

    if (mGuiBeheavior != mBeheavior)
    {
      switch (mGuiBeheavior)
      {
        case 0:
          mNode->getController()->setBeheaviorBalance();
          break;
        case 1:
          mNode->getController()->setBeheaviorWalk();
          break;
        case 2:
          mNode->getController()->setBeheaviorFreeBalance();
          break;
      }

      mBeheavior = mGuiBeheavior;
     }
  }

  if (ImGui::CollapsingHeader("World Options", ImGuiTreeNodeFlags_DefaultOpen))
  {
    // Headlights
    mGuiHeadlights = mViewer->checkHeadlights();
    ImGui::Checkbox("Headlights On/Off", &mGuiHeadlights);
    mViewer->switchHeadlights(mGuiHeadlights);
  }

  if (ImGui::CollapsingHeader("MPC Controller Options", ImGuiTreeNodeFlags_DefaultOpen)) 
  {
    const auto reset = ImGui::Button("Reset Robot");
    if (reset)
      mNode->reset();

    ImGui::Spacing();

    // Gains
    ImGui::SliderFloat( "Vel Gain", &mGuiVelGain, 1, 10000, "%.1g", 10 );
    mNode->getController()->setVelGain(mGuiVelGain);

    ImGui::Spacing();

    ImGui::SliderFloat( "Zmp Pos Gain", &mGuiZmpGain, 1, 10000, "%.1g", 10 );
    mNode->getController()->setZmpGain(mGuiZmpGain);

    ImGui::Spacing();

    // CoM height
    ImGui::SliderFloat("CoM height", &mGuiComHeight, 0.28, 0.33, "%.2f");
    setComTargetHeight(mGuiComHeight);

    ImGui::Spacing();

    // Reference Velocity X
    ImGui::SliderFloat("vref X", &mGuiReferenceVelocityX, -0.2, 0.2, "%.2f");
    setReferenceVelocityX(mGuiReferenceVelocityX);

    ImGui::Spacing();

    // Reference Velocity Y
    ImGui::SliderFloat("vref Y", &mGuiReferenceVelocityY, -0.2, 0.2, "%.2f");
    setReferenceVelocityY(mGuiReferenceVelocityY);

    ImGui::Spacing();

    // Reference Velocity Omega
    ImGui::SliderFloat("vref Omega", &mGuiReferenceVelocityOmega, -0.3, 0.3, "%.2f");
    setReferenceVelocityOmega(mGuiReferenceVelocityOmega);

    ImGui::Spacing();

    // Balance Point
    ImGui::RadioButton("Balance Point: Torso", &mGuiControlMode, 0);
    ImGui::RadioButton("Balance Point: CoM", &mGuiControlMode, 1);

    if (mGuiControlMode != mControlMode)
    {
      switch (mGuiControlMode)
      {
        case 0:
          mNode->getController()->setBalancePointTorso();
          break;
        case 1:
          mNode->getController()->setBalancePointCom();
          break;
        }

      mControlMode = mGuiControlMode;
    }
  }

  if (ImGui::CollapsingHeader("Balance Controller Options", ImGuiTreeNodeFlags_DefaultOpen)) 
  {
    Eigen::VectorXd temp;

    // CoM roll
    ImGui::SliderFloat("CoM roll", &mGuiComRoll, -3.14/2, 3.14/2, "%.3f");
    temp = mNode->getController()->getBalanceFootPos();
    temp(0) = mGuiComRoll;
    mNode->getController()->setBalanceFootPos(temp);

    ImGui::Spacing();

    // CoM pitch
    ImGui::SliderFloat("CoM pitch", &mGuiComPitch, -3.14/2, 3.14/2, "%.3f");
    temp = mNode->getController()->getBalanceFootPos();
    temp(1) = mGuiComPitch;
    mNode->getController()->setBalanceFootPos(temp);

    ImGui::Spacing();

    // CoM yaw
    ImGui::SliderFloat("CoM yaw", &mGuiComYaw, -3.14/2, 3.14/2, "%.3f");
    temp = mNode->getController()->getBalanceFootPos();
    temp(2) = mGuiComYaw;
    mNode->getController()->setBalanceFootPos(temp);

    ImGui::Spacing();

    // CoM x
    ImGui::SliderFloat("CoM X", &mGuiComX, -0.1, 0.1, "%.2f");
    temp = mNode->getController()->getBalanceFootPos();
    temp(3) = mGuiComX;
    mNode->getController()->setBalanceFootPos(temp);

    ImGui::Spacing();

    // CoM y
    ImGui::SliderFloat("CoM Y", &mGuiComY, -0.1, 0.1, "%.2f");
    temp = mNode->getController()->getBalanceFootPos();
    temp(4) = mGuiComY;
    mNode->getController()->setBalanceFootPos(temp);

    ImGui::Spacing();

    // CoM z
    ImGui::SliderFloat("CoM Z", &mGuiComZ, 0.2, 0.4, "%.2f");
    temp = mNode->getController()->getBalanceFootPos();
    temp(5) = mGuiComZ;
    mNode->getController()->setBalanceFootPos(temp);

    ImGui::Spacing();
  }

  // -------------------------- External Forces ---------------------

  if (ImGui::CollapsingHeader("External Forces Options", ImGuiTreeNodeFlags_DefaultOpen)) {
    // Balance Point
    ImGui::RadioButton("External force mode: Constant", &mGuiExternalForceMode, 0);
    ImGui::RadioButton("External force mode: Periodic", &mGuiExternalForceMode, 1);

    if (mGuiExternalForceMode != mExternalForceMode) {
      switch (mGuiBeheavior) {
        case 0:
          mNode->getController()->setExternalForceConstant();
          break;
        case 1:
          mNode->getController()->setExternalForcePeriodic();
          break;
      }
      mExternalForceMode = mGuiExternalForceMode;
    }

    // ImGui::Spacing();
    // // Reference Velocity Omega
    // ImGui::SliderInt("start frame", &mGuiExternalForceStartFrame, 0, 1000, "%d");
    // setExternalForceStartFrame(mGuiExternalForceStartFrame);

    ImGui::Spacing();

    // External force start frame
    ImGui::SliderInt("start frame", &mGuiExternalForceStartFrame, 0, 1000);
    setExternalForceStartFrame(mGuiExternalForceStartFrame);
    
    ImGui::Spacing();

    // External force x
    ImGui::SliderFloat("ext force x", &mGuiExternalForceX, -200.0, 200.0, "%.2f");
    setExternalForceX(mGuiExternalForceX);
    
    ImGui::Spacing();

    // External force y
    ImGui::SliderFloat("ext force y", &mGuiExternalForceY, -200.0, 200.0, "%.2f");
    setExternalForceY(mGuiExternalForceY);
    
    ImGui::Spacing();

    // External force z
    ImGui::SliderFloat("ext force z", &mGuiExternalForceZ, -200.0, 200.0, "%.2f");
    setExternalForceZ(mGuiExternalForceZ);

    ImGui::Spacing();

    // Periodic phase
    ImGui::SliderFloat("extf periodic phase", &mGuiExternalForcePeriodicPhase,
     0.0, 10.0, "%.2f");
    setExternalForcePeriodicPhase(mGuiExternalForcePeriodicPhase);

    ImGui::Spacing();

    // Periodic frequency
    ImGui::SliderFloat("extf periodic freq", &mGuiExternalForcePeriodicFrequency,
     0.0, 2.0, "%.2f");
    setExternalForcePeriodicFrequency(mGuiExternalForcePeriodicFrequency);

    ImGui::Spacing();
  }

  ImGui::End();
}

//==============================================================================
void NaoWidget::setComTargetHeight(double height)
{
  if (mComHeight == height)
    return;

  mComHeight = height;
  mNode->getController()->setComTargetHeight(height);
}

void NaoWidget::setReferenceVelocityX(double ref) {
	mNode->getController()->setReferenceVelocityX(ref);
}

void NaoWidget::setReferenceVelocityY(double ref) {
	mNode->getController()->setReferenceVelocityY(ref);
}

void NaoWidget::setReferenceVelocityOmega(double ref) {
	mNode->getController()->setReferenceVelocityOmega(ref);
}

void NaoWidget::setExternalForceStartFrame(int start)
{
  mNode->getController()->setExternalForceStartFrame(start);
}

void NaoWidget::setExternalForceX(float f)
{
  mNode->getController()->setExternalForceX(f);
}

void NaoWidget::setExternalForceY(float f)
{
  mNode->getController()->setExternalForceY(f);
}

void NaoWidget::setExternalForceZ(float f)
{
  mNode->getController()->setExternalForceZ(f);
}

void NaoWidget::setExternalForcePeriodicPhase(float phi)
{
  mNode->getController()->setExternalForcePeriodicPhase(phi);
}

void NaoWidget::setExternalForcePeriodicFrequency(float k)
{
  mNode->getController()->setExternalForcePeriodicFrequency(k);
}